#' Continuous shrinkage Q95 calibration
#' @export
rc_q95_shrink <- function(C_raw, unit_meta = NULL, stratum_col = NULL, q = 0.95, n0 = 80, eps = 1e-6) {
  C_raw <- as.matrix(C_raw)
  if (is.null(rownames(C_raw)) || is.null(colnames(C_raw))) stop("`C_raw` must have reaction rownames and unit colnames.", call. = FALSE)
  global_q <- apply(C_raw, 1, rc_safe_quantile, probs = q)
  global_n <- rowSums(is.finite(C_raw))
  if (!is.null(stratum_col)) {
    if (is.null(unit_meta) || !"pool_id" %in% colnames(unit_meta) || !stratum_col %in% colnames(unit_meta)) stop("`unit_meta` with `pool_id` and `stratum_col` is required.", call. = FALSE)
    unit_meta <- unit_meta[match(colnames(C_raw), unit_meta$pool_id), , drop = FALSE]
    if (anyNA(unit_meta$pool_id)) stop("`unit_meta` is missing metadata for some capacity columns.", call. = FALSE)
    strata <- as.character(unit_meta[[stratum_col]])
    if (length(unique(stats::na.omit(strata))) < 2L) {
      warning("Only one stratum detected; Q95 stratified calibration degenerates to global calibration.", call. = FALSE)
    }
  } else {
    strata <- rep("global", ncol(C_raw))
  }
  diag_list <- list(); C_rel <- C_raw; idx <- 1L
  for (st in unique(strata)) {
    pools <- which(strata == st)
    n <- rowSums(is.finite(C_raw[, pools, drop = FALSE]))
    rho <- n / (n + n0)
    q_st <- apply(C_raw[, pools, drop = FALSE], 1, rc_safe_quantile, probs = q)
    q_base <- ifelse(is.finite(q_st), q_st, global_q)
    q_shrink <- rho * q_base + (1 - rho) * global_q
    C_rel[, pools] <- sweep(C_raw[, pools, drop = FALSE], 1, q_shrink + eps, "/")
    diag_list[[idx]] <- data.frame(reaction_id = rownames(C_raw), stratum = st,
      n = as.integer(n), n_global = as.integer(global_n), q_stratum = as.numeric(q_st), q_stratum_used = as.numeric(q_base),
      q_global = as.numeric(global_q), rho_n = as.numeric(rho), q_shrink = as.numeric(q_shrink),
      q95_power_class = factor(
        ifelse(n < 5L, "very_low",
          ifelse(n < 20L, "low",
            ifelse(n < 100L, "moderate",
              ifelse(n < 400L, "adequate", "high")
            )
          )
        ),
        levels = c("very_low", "low", "moderate", "adequate", "high"),
        ordered = TRUE
      ),
      stringsAsFactors = FALSE)
    idx <- idx + 1L
  }
  C_rel[C_rel > 1] <- 1
  all_missing <- rowSums(is.finite(C_raw)) == 0L
  if (any(all_missing)) C_rel[all_missing, ] <- NA_real_
  Q <- do.call(rbind, diag_list)
  Q$all_missing_reaction_flag <- Q$n == 0L
  list(C_rel = C_rel, Q = Q)
}

#' Calibrate raw reaction capacities by continuous reaction-wise Q95 shrinkage
#'
#' Continuous shrinkage is used for every reaction/stratum with one or more
#' finite units; low unit counts are reported through diagnostics rather than a
#' hard direct/global switch.
#' @export
rc_q95_calibrate <- function(C_raw, eps = 1e-6, bootstrap = TRUE, B = 500, BPPARAM = NULL, n0 = 80, unit_meta = NULL, stratum_col = NULL) {
  C_raw <- as.matrix(C_raw)
  out <- rc_q95_shrink(C_raw, unit_meta = unit_meta, stratum_col = stratum_col, q = 0.95, n0 = n0, eps = eps)
  Q <- out$Q
  names(Q)[names(Q) == "q_shrink"] <- "q_value"
  Q$quantile_used <- 0.95
  Q$n_finite <- Q$n
  Q$n_finite_global <- rowSums(is.finite(C_raw))[match(Q$reaction_id, rownames(C_raw))]
  Q$low_n_flag <- Q$n_finite < 20L
  if (isTRUE(bootstrap)) {
    boot <- rc_q95_bootstrap_diagnostics(C_raw, Q, unit_meta = unit_meta, stratum_col = stratum_col, B = B, BPPARAM = BPPARAM)
    Q$q95_bootstrap <- boot[, "q95"]; Q$q95_ci_low <- boot[, "ci_low"]
    Q$q95_ci_high <- boot[, "ci_high"]; Q$q95_ci_width <- boot[, "width"]
    Q$q95_unstable_flag <- is.finite(Q$q95_ci_width) &
      (Q$q95_ci_width / pmax(Q$q_value, 1e-6)) > 0.5
  }
  list(C_rel = out$C_rel, Q = Q)
}

#' Calibrate raw reaction capacity by Q95 shrinkage
#' @export
rc_q95_bootstrap <- function(x, B = 500) {
  x <- x[is.finite(x)]
  if (length(x) < 20L) {
    return(c(q95 = NA_real_, ci_low = NA_real_, ci_high = NA_real_, width = NA_real_))
  }
  if (!is.numeric(B) || length(B) != 1L || is.na(B) || B < 1) {
    stop("`B` must be a single positive number.", call. = FALSE)
  }
  qs <- replicate(B, stats::quantile(sample(x, replace = TRUE), 0.95, na.rm = TRUE, names = FALSE))
  ci <- stats::quantile(qs, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  c(
    q95 = stats::quantile(x, 0.95, na.rm = TRUE, names = FALSE),
    ci_low = unname(ci[[1]]),
    ci_high = unname(ci[[2]]),
    width = unname(diff(ci))
  )
}


rc_empty_gpr_evidence_matrix <- function(gpr_list, unit_ids) {
  genes <- unique(tolower(unlist(gpr_list, use.names = FALSE)))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (length(genes) == 0L) genes <- character(0)
  matrix(NA_real_, nrow = length(genes), ncol = length(unit_ids), dimnames = list(genes, unit_ids))
}

rc_set_confidence_source <- function(confidence, source, detection_available = NULL) {
  confidence$confidence_source <- source
  confidence$evidence_source <- source
  if (!is.null(detection_available)) confidence$detection_available <- detection_available
  confidence
}

#' Run the simplified Layer 1 capacity workflow
#' @export
rc_run_layer1_capacity <- function(gpr_table,
                                   unit_expression,
                                   unit_detection = NULL,
                                   unit_meta = NULL,
                                   stratum_col = NULL,
                                   gene_confidence = NULL,
                                   gene_discordance = NULL,
                                   and_methods = c("min", "boltzmann_0.08", "boltzmann_0.20", "mean"),
                                   tau_sensitivity = c(0.08, 0.20),
                                   promiscuity_sensitivity = c("none", "sqrt", "linear"),
                                   promiscuity_mode = c("sqrt", "linear", "none"),
                                   tau = 0.20,
                                   and_method = c("boltzmann", "min", "mean"),
                                   or_method = c("sum_sqrtK", "max", "prob_or", "sum"),
                                   run_sensitivity = FALSE,
                                   bootstrap = FALSE,
                                   low_confidence_threshold = 0.25,
                                   low_confidence_quantile = NULL,
                                   reaction_confidence_method = c("gpr_aware"),
                                   B = 500,
                                   BPPARAM = NULL) {
  promiscuity_mode <- match.arg(promiscuity_mode)
  and_method <- match.arg(and_method)
  or_method <- match.arg(or_method)
  reaction_confidence_method <- match.arg(reaction_confidence_method)
  unit_expression <- as.matrix(unit_expression)
  if (is.null(rownames(unit_expression))) {
    stop("`unit_expression` must have gene IDs in rownames().", call. = FALSE)
  }
  if (is.null(colnames(unit_expression))) {
    stop("`unit_expression` must have unit IDs in colnames().", call. = FALSE)
  }
  parsed <- rc_parse_gpr_table(gpr_table)
  gene_score <- rc_gene_score(unit_expression)
  rownames(gene_score) <- tolower(rownames(gene_score))
  if (!is.null(unit_detection)) {
    unit_detection <- as.matrix(unit_detection)
    if (is.null(rownames(unit_detection)) || is.null(colnames(unit_detection))) stop("`unit_detection` must have gene rownames and unit colnames.", call. = FALSE)
    if (!all(colnames(unit_expression) %in% colnames(unit_detection))) stop("`unit_detection` is missing one or more unit_expression columns.", call. = FALSE)
    unit_detection <- unit_detection[, colnames(unit_expression), drop = FALSE]
  }
  if (!is.null(gene_confidence)) {
    gene_confidence <- as.matrix(gene_confidence)
    if (is.null(rownames(gene_confidence)) || is.null(colnames(gene_confidence))) stop("`gene_confidence` must have gene rownames and unit colnames.", call. = FALSE)
    if (!all(colnames(unit_expression) %in% colnames(gene_confidence))) stop("`gene_confidence` is missing one or more unit_expression columns.", call. = FALSE)
    gene_confidence <- gene_confidence[, colnames(unit_expression), drop = FALSE]
  }

  C_raw <- rc_reaction_capacity(parsed, gene_score, promiscuity_mode = promiscuity_mode, tau = tau, and_method = and_method, or_method = or_method, BPPARAM = BPPARAM)
  calibrated <- rc_q95_calibrate(C_raw, bootstrap = bootstrap, B = B, BPPARAM = BPPARAM, unit_meta = unit_meta, stratum_col = stratum_col)
  gpr_diag <- rc_gpr_diagnostics(parsed, rownames(gene_score))
  confidence <- rc_reaction_confidence(
    parsed,
    gene_confidence = gene_confidence,
    unit_detection = unit_detection,
    method = reaction_confidence_method,
    tau_conf = tau,
    low_confidence_quantile = low_confidence_quantile,
    low_confidence_threshold = low_confidence_threshold,
    unit_ids = colnames(unit_expression)
  )
  confidence$reaction_confidence_method <- reaction_confidence_method
  if (isTRUE(run_sensitivity)) {
    tau_sens <- rc_capacity_sensitivity(parsed, gene_score, variable = "tau", values = tau_sensitivity, promiscuity_mode = promiscuity_mode, and_method = and_method, or_method = or_method, BPPARAM = BPPARAM)
    prom_sens <- rc_capacity_sensitivity(parsed, gene_score, variable = "promiscuity", values = promiscuity_sensitivity, tau = tau, and_method = and_method, or_method = or_method, BPPARAM = BPPARAM)
    and_long <- rc_and_method_capacity_long(parsed, gene_score, and_methods = and_methods, promiscuity_mode = promiscuity_mode, or_method = or_method, BPPARAM = BPPARAM)
    and_sens <- rc_and_method_sensitivity(and_long)
  } else {
    tau_sens <- data.frame()
    prom_sens <- data.frame()
    and_long <- data.frame()
    and_sens <- data.frame()
  }

  confidence_sources <- rc_reaction_confidence_sources(confidence)
  confidence_summary <- rc_reaction_confidence_summary(confidence, gene_confidence = gene_confidence)

  out <- list(
    C_or_raw = C_raw,
    or_method_used = or_method,
    C_raw = C_raw,
    reaction_capacity_L1 = C_raw,
    C_rel = calibrated$C_rel,
    reaction_confidence = confidence,
    reaction_confidence_source = confidence_sources,
    reaction_confidence_summary = confidence_summary,
    capacity_long = and_long,
    and_method_sensitivity = and_sens,
    q95_diagnostics = calibrated$Q,
    gpr_diagnostics = gpr_diag,
    gene_discordance = gene_discordance,
    tau_sensitivity = tau_sens,
    promiscuity_sensitivity = prom_sens,
    parsed_gpr = parsed,
    reaction_confidence_method = reaction_confidence_method,
    unit_meta = unit_meta
  )
  if (identical(or_method, "sum")) out <- c(list(C_iso_sum_raw = C_raw), out)
  out
}


#' Soft-min aggregation for confidence scores
rc_softmin_conf <- function(v, tau = 0.20) {
  v <- as.numeric(v)
  if (any(!is.finite(v))) return(NA_real_)
  w <- exp(-v / tau)
  sum(w * v) / sum(w)
}

rc_apply_low_confidence_quantile <- function(df, low_confidence_quantile) {
  df$low_confidence_reaction_flag <- NA
  if (is.null(low_confidence_quantile)) return(df)
  if (!is.numeric(low_confidence_quantile) || length(low_confidence_quantile) != 1L || is.na(low_confidence_quantile) || low_confidence_quantile <= 0 || low_confidence_quantile >= 1) {
    stop("`low_confidence_quantile` must be NULL or a single number between 0 and 1.", call. = FALSE)
  }
  split_src <- split(seq_len(nrow(df)), df$confidence_source)
  for (idx in split_src) {
    x <- df$reaction_confidence[idx]
    if (!any(is.finite(x))) next
    threshold <- stats::quantile(x, probs = low_confidence_quantile, na.rm = TRUE, names = FALSE)
    df$low_confidence_reaction_flag[idx] <- is.finite(x) & x < threshold
  }
  df
}

#' Compute GPR-aware reaction-level confidence from gene confidence
#'
#' AND groups are aggregated with softmin/min/mean, OR groups with max/prob_or/sum_sqrtK.
#' Missing required genes make the affected AND group incomplete; a reaction remains
#' supported when at least one alternative AND group is complete.
#' @export
rc_reaction_confidence_gpr_aware <- function(gpr_list,
                                             gene_confidence,
                                             unit_detection = NULL,
                                             tau_conf = 0.20,
                                             and_method = c("softmin", "min", "mean"),
                                             or_method = c("max", "prob_or", "sum_sqrtK"),
                                             missing_group_policy = c("complete_group", "partial_group"),
                                             low_confidence_quantile = NULL) {
  and_method <- match.arg(and_method)
  or_method <- match.arg(or_method)
  missing_group_policy <- match.arg(missing_group_policy)
  gene_confidence <- as.matrix(gene_confidence)
  if (is.null(rownames(gene_confidence)) || is.null(colnames(gene_confidence))) stop("`gene_confidence` must have gene rownames and unit colnames.", call. = FALSE)
  rownames(gene_confidence) <- tolower(rownames(gene_confidence))
  unit_ids <- colnames(gene_confidence)

  rows <- lapply(names(gpr_list), function(rid) {
    groups <- gpr_list[[rid]]
    group_mat <- lapply(seq_along(groups), function(k) {
      gs <- unique(tolower(groups[[k]]))
      gs <- gs[!is.na(gs) & nzchar(gs)]
      hit <- intersect(gs, rownames(gene_confidence))
      observed_fraction <- if (length(gs) == 0L) NA_real_ else length(hit) / length(gs)
      complete <- length(gs) > 0L && length(hit) == length(gs)
      if (!complete && identical(missing_group_policy, "complete_group")) {
        return(list(conf = rep(NA_real_, length(unit_ids)), observed_fraction = observed_fraction, complete = FALSE))
      }
      if (length(hit) == 0L) {
        return(list(conf = rep(NA_real_, length(unit_ids)), observed_fraction = observed_fraction, complete = FALSE))
      }
      x <- gene_confidence[hit, , drop = FALSE]
      conf <- if (nrow(x) == 1L) {
        as.numeric(x[1, ])
      } else if (identical(and_method, "min")) {
        apply(x, 2, function(v) if (any(!is.finite(v))) NA_real_ else min(v))
      } else if (identical(and_method, "mean")) {
        colMeans(x, na.rm = FALSE)
      } else {
        apply(x, 2, rc_softmin_conf, tau = tau_conf)
      }
      list(conf = conf, observed_fraction = observed_fraction, complete = complete)
    })
    if (length(group_mat) == 0L) group_mat <- list(list(conf = rep(NA_real_, length(unit_ids)), observed_fraction = NA_real_, complete = FALSE))
    G <- do.call(rbind, lapply(group_mat, `[[`, "conf"))
    complete <- vapply(group_mat, `[[`, logical(1), "complete")
    observed <- vapply(group_mat, `[[`, numeric(1), "observed_fraction")
    if (!any(complete)) {
      vals <- rep(NA_real_, length(unit_ids))
    } else {
      G2 <- G[complete, , drop = FALSE]
      vals <- if (identical(or_method, "max")) {
        apply(G2, 2, function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE))
      } else if (identical(or_method, "prob_or")) {
        apply(G2, 2, function(v) { v <- v[is.finite(v)]; if (!length(v)) return(NA_real_); 1 - prod(1 - pmax(0, pmin(1, v))) })
      } else {
        apply(G2, 2, function(v) { v <- v[is.finite(v)]; if (!length(v)) return(NA_real_); sum(v) / sqrt(length(v)) })
      }
    }
    genes_all <- unique(tolower(unlist(groups, use.names = FALSE)))
    genes_all <- genes_all[!is.na(genes_all) & nzchar(genes_all)]
    data.frame(
      reaction_id = rid, pool_id = unit_ids,
      reaction_confidence = vals, reaction_evidence_score = vals,
      reaction_confidence_unpenalized = vals, reaction_evidence_score_unpenalized = vals,
      confidence_source = "gpr_aware_gene_confidence", evidence_source = "gpr_aware_gene_confidence",
      n_gpr_genes_total = length(genes_all), n_gpr_genes_multiome = length(intersect(genes_all, rownames(gene_confidence))),
      multiome_coverage_fraction = if (length(genes_all) == 0L) NA_real_ else length(intersect(genes_all, rownames(gene_confidence))) / length(genes_all),
      missing_gpr_gene_fraction = if (length(genes_all) == 0L) NA_real_ else 1 - length(intersect(genes_all, rownames(gene_confidence))) / length(genes_all),
      missing_gene_fraction = if (length(genes_all) == 0L) NA_real_ else 1 - length(intersect(genes_all, rownames(gene_confidence))) / length(genes_all),
      missing_subunit_confidence_penalty = NA_real_, detection_available = !is.null(unit_detection), mean_gpr_detection_rate = NA_real_,
      low_confidence_threshold = NA_real_,
      n_and_groups_total = length(groups), n_and_groups_complete = sum(complete),
      complete_and_group_fraction = mean(complete),
      best_and_group_observed_fraction = if (all(is.na(observed))) NA_real_ else max(observed, na.rm = TRUE),
      any_incomplete_gpr_group_flag = any(!complete),
      reaction_unsupported_by_complete_gpr_flag = !any(complete),
      missing_required_subunit_flag = any(!complete),
      no_complete_gpr_group_flag = !any(complete),
      stringsAsFactors = FALSE
    )
  })
  rc_apply_low_confidence_quantile(do.call(rbind, rows), low_confidence_quantile)
}


#' Aggregate one evidence matrix through GPR AND/OR structure
#'
#' Confidence is bottleneck-aware for AND groups (minimum subunit evidence) and
#' isoenzyme-aware for OR groups (maximum supported isoenzyme evidence). Missing
#' genes are handled by the reaction-level missing-gene penalty rather than being silently ignored.
rc_gpr_confidence_one <- function(parsed_gpr, evidence_vec) {
  and_vals <- vapply(parsed_gpr, function(and_group) {
    genes <- unique(tolower(and_group))
    vals <- evidence_vec[genes]
    vals <- vals[is.finite(vals)]
    if (length(genes) == 0L || length(vals) == 0L) {
      return(NA_real_)
    }
    min(vals)
  }, numeric(1))
  finite <- and_vals[is.finite(and_vals)]
  if (length(finite) == 0L) return(NA_real_)
  max(finite)
}

rc_gpr_confidence_matrix <- function(gpr, evidence) {
  vapply(seq_len(ncol(evidence)), function(j) {
    rc_gpr_confidence_one(gpr, evidence[, j])
  }, numeric(1))
}

#' Compute reaction-level confidence from gene confidence or RNA detection
#'
#' By default this public API uses GPR-aware aggregation for either multiome
#' gene confidence or RNA-only detection evidence.
#' @export
rc_reaction_confidence <- function(gpr_list,
                                   gene_confidence = NULL,
                                   unit_detection = NULL,
                                   method = c("gpr_aware"),
                                   tau_conf = 0.20,
                                   and_method = c("softmin", "min", "mean"),
                                   or_method = c("max", "prob_or", "sum_sqrtK"),
                                   low_confidence_quantile = NULL,
                                   low_confidence_threshold = NULL,
                                   unit_ids = NULL) {
  method <- match.arg(method)
  and_method <- match.arg(and_method)
  or_method <- match.arg(or_method)

  if (identical(method, "gpr_aware")) {
    evidence <- if (!is.null(gene_confidence)) gene_confidence else unit_detection
    source <- if (!is.null(gene_confidence)) "gpr_aware_gene_confidence" else "gpr_aware_rna_detection"
    detection_available <- !is.null(unit_detection)
    if (is.null(evidence)) {
      if (is.null(unit_ids)) stop("Provide `gene_confidence` or `unit_detection` for gpr_aware confidence.", call. = FALSE)
      evidence <- rc_empty_gpr_evidence_matrix(gpr_list, unit_ids)
      source <- "gpr_aware_no_evidence"
      detection_available <- FALSE
    }
    out <- rc_reaction_confidence_gpr_aware(
      gpr_list = gpr_list,
      gene_confidence = evidence,
      unit_detection = unit_detection,
      tau_conf = tau_conf,
      and_method = and_method,
      or_method = or_method,
      low_confidence_quantile = low_confidence_quantile
    )
    return(rc_set_confidence_source(out, source, detection_available = detection_available))
  }
}

rc_gpr_gene_ids <- function(gpr_list) {
  genes <- unique(tolower(unlist(gpr_list, use.names = FALSE)))
  genes[!is.na(genes) & nzchar(genes)]
}

#' Extract metabolic GPR genes from a GPR table
#'
#' Use this gene set as the `genes.use` target when regenerating metabolic
#' peak-gene links with Signac::LinkPeaks for multiome confidence.
#' @export
rc_metabolic_gpr_genes <- function(gpr_table) {
  if (is.data.frame(gpr_table) && all(c("reaction_id", "and_group_id", "gene") %in% colnames(gpr_table))) {
    genes <- unique(as.character(gpr_table$gene))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    return(genes)
  }
  if (is.data.frame(gpr_table) && all(c("reaction_id", "gpr") %in% colnames(gpr_table))) {
    rules <- as.character(gpr_table$gpr)
    rules <- gsub("[()]", " ", rules)
    rules <- gsub("\\s+(and|or)\\s+", "|", rules, perl = TRUE, ignore.case = TRUE)
    tokens <- unlist(strsplit(rules, "|", fixed = TRUE), use.names = FALSE)
    genes <- unique(trimws(tokens))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    return(genes)
  }
  gpr_list <- if (is.data.frame(gpr_table)) rc_parse_gpr_table(gpr_table) else gpr_table
  genes <- unique(as.character(unlist(gpr_list, use.names = FALSE)))
  genes[!is.na(genes) & nzchar(genes)]
}

#' Compute compact capacity sensitivity summaries
#' @export
rc_capacity_sensitivity <- function(gpr_list, gene_score, variable = c("tau", "promiscuity"), values, promiscuity_mode = "sqrt", tau = 0.20, and_method = "boltzmann", or_method = "sum", BPPARAM = NULL) {
  variable <- match.arg(variable)
  if (length(values) == 0L) return(data.frame())
  rows <- lapply(values, function(v) {
    C <- if (variable == "tau") {
      rc_reaction_capacity(gpr_list, gene_score, promiscuity_mode = promiscuity_mode, tau = as.numeric(v), and_method = and_method, or_method = or_method, BPPARAM = BPPARAM)
    } else {
      rc_reaction_capacity(gpr_list, gene_score, promiscuity_mode = as.character(v), tau = tau, and_method = and_method, or_method = or_method, BPPARAM = BPPARAM)
    }
    data.frame(reaction_id = rownames(C), sensitivity = variable, value = as.character(v), median_capacity = matrixStats::rowMedians(C, na.rm = TRUE), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

rc_safe_quantile <- function(x, probs) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
}


rc_q95_bootstrap_diagnostics <- function(C_raw, Q, unit_meta = NULL, stratum_col = NULL, B = 500, BPPARAM = NULL) {
  C_raw <- as.matrix(C_raw)
  strata <- if (!is.null(stratum_col)) {
    unit_meta <- unit_meta[match(colnames(C_raw), unit_meta$pool_id), , drop = FALSE]
    as.character(unit_meta[[stratum_col]])
  } else rep("global", ncol(C_raw))
  rows <- seq_len(nrow(Q))
  boot <- rc_parallel_lapply(rows, function(i) {
    rid <- Q$reaction_id[[i]]; st <- Q$stratum[[i]]
    cols <- which(strata == st)
    rc_q95_bootstrap(C_raw[rid, cols], B = B)
  }, BPPARAM = BPPARAM)
  do.call(rbind, boot)
}

#' Layer 1 capacity alias matching the adjusted plan examples
#' @export
rc_layer1_capacity <- rc_run_layer1_capacity

#' Compute long-form capacity for multiple AND methods
#' @export
rc_and_method_capacity_long <- function(gpr_list, gene_score, and_methods = c("min", "boltzmann_0.08", "boltzmann_0.20", "mean"), promiscuity_mode = "sqrt", or_method = "sum", BPPARAM = NULL) {
  rows <- lapply(and_methods, function(m) {
    spec <- rc_parse_and_method(m)
    C <- rc_reaction_capacity(gpr_list, gene_score, promiscuity_mode = promiscuity_mode, tau = spec$tau, and_method = spec$method, or_method = or_method, BPPARAM = BPPARAM)
    data.frame(reaction_id = rep(rownames(C), times = ncol(C)), pool_id = rep(colnames(C), each = nrow(C)), and_method = m, C_raw = as.vector(C), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Summarize AND-method sensitivity by reaction and pool
#' @export
rc_and_method_sensitivity <- function(capacity_long) {
  if (nrow(capacity_long) == 0L) return(data.frame())
  split_key <- interaction(capacity_long$reaction_id, capacity_long$pool_id, drop = TRUE)
  rows <- lapply(split(capacity_long, split_key), function(df) {
    tau_vals <- df$C_raw[df$and_method %in% c("boltzmann_0.08", "boltzmann_0.20")]
    tau_vals <- tau_vals[is.finite(tau_vals)]
    data.frame(reaction_id = df$reaction_id[[1]], pool_id = df$pool_id[[1]],
               capacity_range = rc_safe_range(df$C_raw),
               tau_sensitive_flag = length(tau_vals) >= 2L && diff(range(tau_vals)) > 0.05,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

rc_safe_range <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  max(x) - min(x)
}

rc_parse_and_method <- function(x) {
  if (x == "min") return(list(method = "min", tau = 0.20))
  if (x == "mean") return(list(method = "mean", tau = 0.20))
  if (grepl("^boltzmann_", x)) return(list(method = "boltzmann", tau = as.numeric(sub("^boltzmann_", "", x))))
  stop("Unsupported AND method: ", x, call. = FALSE)
}

#' End-to-end Layer 1 run from raw counts and a pool map
rc_run_layer1_from_counts <- function(gpr_table,
                                      rna_counts,
                                      unit_map,
                                      unit_meta = NULL,
                                      atac_counts = NULL,
                                      peak_gene_links = NULL,
                                      stratum_col = "cell_type",
                                      promiscuity_mode = "sqrt",
                                      and_method = "boltzmann",
                                      tau = 0.20,
                                      or_method = c("sum_sqrtK", "max", "prob_or", "sum"),
                                      bootstrap = FALSE,
                                      low_confidence_threshold = 0.25,
                                      low_confidence_quantile = NULL,
                                      reaction_confidence_method = c("gpr_aware"),
                                      B = 500,
                                      BPPARAM = NULL) {
  or_method <- match.arg(or_method)
  reaction_confidence_method <- match.arg(reaction_confidence_method)
  rna_pb <- rc_unit_bulk_counts(rna_counts, unit_map, fun = "sum", BPPARAM = BPPARAM)
  if (is.null(unit_meta)) unit_meta <- rc_build_unit_metadata(unit_map)
  filtered <- rc_filter_empty_units(rna_pb, unit_meta)
  rna_logcpm <- rc_logcpm(filtered$counts)
  unit_meta <- filtered$unit_meta
  rna_detection <- rc_unit_detection(rna_counts, unit_map, BPPARAM = BPPARAM)
  rna_detection <- rna_detection[, colnames(rna_logcpm), drop = FALSE]

  gene_conf <- NULL
  gene_confidence_components <- NULL
  parsed_gpr <- rc_parse_gpr_table(gpr_table)
  metabolic_gpr_genes <- rc_metabolic_gpr_genes(parsed_gpr)
  if (!is.null(atac_counts) && !is.null(peak_gene_links)) {
    peak_gene_links <- rc_filter_peak_gene_links_to_gpr(peak_gene_links, metabolic_gpr_genes)
    atac_peak <- rc_atac_unit_logcpm(atac_counts, unit_map, min_pools = 3, BPPARAM = BPPARAM)
    common_pools <- intersect(colnames(rna_logcpm), colnames(atac_peak))
    rna_logcpm <- rna_logcpm[, common_pools, drop = FALSE]
    rna_detection <- rna_detection[, common_pools, drop = FALSE]
    unit_meta <- unit_meta[match(common_pools, unit_meta$pool_id), , drop = FALSE]
    p_rna <- rc_percentile_by_stratum(rna_logcpm, unit_meta = unit_meta, stratum_col = stratum_col)
    p_atac_peak <- rc_percentile_by_stratum(atac_peak[, common_pools, drop = FALSE], unit_meta = unit_meta, stratum_col = stratum_col)
    if (nrow(peak_gene_links) > 0L) {
      link_conf <- rc_link_confidence(p_atac_peak, peak_gene_links)
      genes <- Reduce(intersect, list(tolower(rownames(p_rna)), tolower(rownames(link_conf)), metabolic_gpr_genes))
    } else {
      genes <- character(0)
    }
    if (length(genes) > 0L) {
      rownames(p_rna) <- tolower(rownames(p_rna))
      rownames(rna_logcpm) <- tolower(rownames(rna_logcpm))
      rownames(rna_detection) <- tolower(rownames(rna_detection))
      rownames(link_conf) <- tolower(rownames(link_conf))
      concord <- rc_concordance_null_correct(p_rna[genes, , drop = FALSE], link_conf[genes, , drop = FALSE], unit_meta = unit_meta, stratum_col = stratum_col)
      rel <- rc_fisher_shrink(rna_logcpm[genes, , drop = FALSE], link_conf[genes, , drop = FALSE])$rel_positive
      names(rel) <- genes
      gene_conf_out <- rc_gene_confidence(concord, rel_ra_pos = rel, det_rna = rna_detection[genes, , drop = FALSE], link_conf = link_conf[genes, , drop = FALSE], return_components = TRUE)
      gene_conf <- gene_conf_out$gene_confidence
      gene_confidence_components <- gene_conf_out[setdiff(names(gene_conf_out), c("gene_confidence", "confidence"))]
    }
  }

  out <- rc_run_layer1_capacity(gpr_table = gpr_table, unit_expression = rna_logcpm, unit_detection = rna_detection,
                                unit_meta = unit_meta, stratum_col = stratum_col, gene_confidence = gene_conf,
                                promiscuity_mode = promiscuity_mode, and_method = and_method, or_method = or_method, tau = tau, bootstrap = bootstrap,
                                low_confidence_threshold = low_confidence_threshold, low_confidence_quantile = low_confidence_quantile,
                                reaction_confidence_method = reaction_confidence_method, B = B, BPPARAM = BPPARAM)
  out$unit_meta <- unit_meta
  out$gene_confidence_components <- gene_confidence_components
  out
}

rc_reaction_confidence_sources <- function(reaction_confidence) {
  if (is.null(reaction_confidence) || nrow(reaction_confidence) == 0L || !"confidence_source" %in% colnames(reaction_confidence)) return(character(0))
  src <- unique(as.character(reaction_confidence$confidence_source))
  src <- src[!is.na(src) & nzchar(src)]
  if (length(src) == 1L && identical(src, "rna_detection_fallback")) return("rna_detection")
  src
}

#' Summarize reaction confidence evidence coverage
#' @export
rc_reaction_confidence_summary <- function(reaction_confidence, gene_confidence = NULL, low_multiome_threshold = 0.20) {
  if (is.null(reaction_confidence) || nrow(reaction_confidence) == 0L) return(data.frame())
  src <- as.character(reaction_confidence$confidence_source)
  rxn_src <- stats::aggregate(src %in% c("multiome_link_confidence", "gpr_aware_gene_confidence"), list(reaction_id = reaction_confidence$reaction_id), any)
  names(rxn_src)[2] <- "has_multiome"
  multiome_reaction_fraction <- mean(rxn_src$has_multiome, na.rm = TRUE)
  multiome_reaction_pool_fraction <- mean(src %in% c("multiome_link_confidence", "gpr_aware_gene_confidence"), na.rm = TRUE)
  fallback_fraction <- mean(src == "rna_detection_fallback", na.rm = TRUE)
  out <- data.frame(
    multiome_coverage_reaction_fraction = multiome_reaction_fraction,
    multiome_coverage_reaction_pool_fraction = multiome_reaction_pool_fraction,
    fallback_fraction = fallback_fraction,
    n_linked_metabolic_genes = if (!is.null(gene_confidence)) nrow(gene_confidence) else NA_integer_,
    stringsAsFactors = FALSE
  )
  if (is.finite(multiome_reaction_pool_fraction) && multiome_reaction_pool_fraction > 0 && multiome_reaction_pool_fraction < low_multiome_threshold) {
    warning("Multiome coverage is low; results are mostly driven by RNA fallback.", call. = FALSE)
  }
  out
}

#' Filter peak-gene links to metabolic GPR genes
#' @export
rc_filter_peak_gene_links_to_gpr <- function(peak_gene_links, gpr_genes) {
  required <- c("peak_id", "gene", "weight")
  missing <- setdiff(required, colnames(peak_gene_links))
  if (length(missing) > 0L) stop("`peak_gene_links` is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  if (length(gpr_genes) == 0L) return(peak_gene_links[0, , drop = FALSE])
  gpr_genes <- unique(tolower(as.character(gpr_genes)))
  peak_gene_links[tolower(as.character(peak_gene_links$gene)) %in% gpr_genes, , drop = FALSE]
}
