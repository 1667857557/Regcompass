rc_q95_shrink <- function(C_raw, unit_meta = NULL, stratum_col = NULL,
                          q = 0.95, n0 = 80, eps = 1e-6) {
  C_raw <- as.matrix(C_raw)
  if (is.null(rownames(C_raw)) || is.null(colnames(C_raw))) {
    stop("`C_raw` must have reaction rownames and unit colnames.", call. = FALSE)
  }
  global_q <- apply(C_raw, 1, rc_safe_quantile, probs = q)
  global_n <- rowSums(is.finite(C_raw))

  if (!is.null(stratum_col)) {
    if (is.null(unit_meta) || !"pool_id" %in% colnames(unit_meta) ||
        !stratum_col %in% colnames(unit_meta)) {
      stop("`unit_meta` with `pool_id` and `stratum_col` is required.",
           call. = FALSE)
    }
    unit_meta <- unit_meta[match(colnames(C_raw), unit_meta$pool_id), , drop = FALSE]
    if (anyNA(unit_meta$pool_id)) {
      stop("`unit_meta` is missing metadata for some capacity columns.",
           call. = FALSE)
    }
    strata <- trimws(as.character(unit_meta[[stratum_col]]))
    if (anyNA(strata) || any(!nzchar(strata))) {
      stop("Q95 strata must be non-missing and non-empty.", call. = FALSE)
    }
    if (length(unique(strata)) < 2L) {
      warning(
        "Only one stratum detected; Q95 stratified calibration degenerates to global calibration.",
        call. = FALSE
      )
    }
  } else {
    strata <- rep("global", ncol(C_raw))
  }

  diagnostics <- list()
  C_rel <- C_raw
  index <- 1L
  for (stratum in unique(strata)) {
    pools <- which(strata == stratum)
    n <- rowSums(is.finite(C_raw[, pools, drop = FALSE]))
    rho <- n / (n + n0)
    q_stratum <- apply(
      C_raw[, pools, drop = FALSE],
      1,
      rc_safe_quantile,
      probs = q
    )
    q_used <- ifelse(is.finite(q_stratum), q_stratum, global_q)
    q_shrink <- rho * q_used + (1 - rho) * global_q
    C_rel[, pools] <- sweep(
      C_raw[, pools, drop = FALSE],
      1,
      q_shrink + eps,
      "/"
    )
    diagnostics[[index]] <- data.frame(
      reaction_id = rownames(C_raw),
      stratum = stratum,
      n = as.integer(n),
      n_global = as.integer(global_n),
      q_stratum = as.numeric(q_stratum),
      q_stratum_used = as.numeric(q_used),
      q_global = as.numeric(global_q),
      rho_n = as.numeric(rho),
      q_shrink = as.numeric(q_shrink),
      q95_power_class = factor(
        ifelse(
          n < 5L,
          "very_low",
          ifelse(
            n < 20L,
            "low",
            ifelse(n < 100L, "moderate",
                   ifelse(n < 400L, "adequate", "high"))
          )
        ),
        levels = c("very_low", "low", "moderate", "adequate", "high"),
        ordered = TRUE
      ),
      stringsAsFactors = FALSE
    )
    index <- index + 1L
  }

  C_rel[C_rel > 1] <- 1
  all_missing <- global_n == 0L
  if (any(all_missing)) C_rel[all_missing, ] <- NA_real_
  Q <- do.call(rbind, diagnostics)
  Q$all_missing_reaction_flag <- Q$n_global == 0L
  Q$stratum_missing_reaction_flag <- Q$n == 0L
  list(C_rel = C_rel, Q = Q)
}

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

rc_reaction_confidence_gpr_aware <- function(
    gpr_list,
    gene_confidence,
    unit_detection = NULL,
    tau_conf = 0.20,
    and_method = c("softmin", "min", "mean"),
    or_method = c("max", "prob_or", "sum_sqrtK"),
    missing_group_policy = c("complete_group", "partial_group"),
    low_confidence_quantile = NULL,
    low_confidence_threshold = NULL) {
  and_method <- match.arg(and_method)
  or_method <- match.arg(or_method)
  missing_group_policy <- match.arg(missing_group_policy)
  if (!is.null(low_confidence_threshold)) {
    if (!is.numeric(low_confidence_threshold) ||
        length(low_confidence_threshold) != 1L ||
        !is.finite(low_confidence_threshold) ||
        low_confidence_threshold < 0 ||
        low_confidence_threshold > 1) {
      stop(
        "`low_confidence_threshold` must be NULL or one number between 0 and 1.",
        call. = FALSE
      )
    }
    if (!is.null(low_confidence_quantile)) {
      stop(
        "Specify only one of `low_confidence_threshold` and `low_confidence_quantile`.",
        call. = FALSE
      )
    }
  }
  gene_confidence <- as.matrix(gene_confidence)
  if (is.null(rownames(gene_confidence)) ||
      is.null(colnames(gene_confidence))) {
    stop(
      "`gene_confidence` must have gene rownames and unit colnames.",
      call. = FALSE
    )
  }
  normalized <- tolower(trimws(rownames(gene_confidence)))
  if (anyNA(normalized) || any(!nzchar(normalized)) ||
      anyDuplicated(normalized)) {
    stop(
      "`gene_confidence` must have unique non-empty genes after normalization.",
      call. = FALSE
    )
  }
  rownames(gene_confidence) <- normalized
  unit_ids <- colnames(gene_confidence)

  rows <- lapply(names(gpr_list), function(rid) {
    groups <- gpr_list[[rid]]
    group_data <- lapply(groups, function(group) {
      genes <- unique(tolower(trimws(as.character(group))))
      genes <- genes[!is.na(genes) & nzchar(genes)]
      hit <- intersect(genes, rownames(gene_confidence))
      observed <- if (!length(genes)) NA_real_ else length(hit) / length(genes)
      complete <- length(genes) > 0L && length(hit) == length(genes)
      eligible <- complete ||
        (identical(missing_group_policy, "partial_group") && length(hit) > 0L)
      if (!eligible) {
        return(list(
          values = rep(NA_real_, length(unit_ids)),
          complete = complete,
          eligible = FALSE,
          observed = observed
        ))
      }
      matrix_in <- gene_confidence[hit, , drop = FALSE]
      values <- if (nrow(matrix_in) == 1L) {
        as.numeric(matrix_in[1L, ])
      } else if (identical(and_method, "min")) {
        apply(matrix_in, 2L, function(x) {
          if (any(!is.finite(x))) NA_real_ else min(x)
        })
      } else if (identical(and_method, "mean")) {
        colMeans(matrix_in, na.rm = FALSE)
      } else {
        apply(matrix_in, 2L, rc_softmin_conf, tau = tau_conf)
      }
      list(
        values = values,
        complete = complete,
        eligible = TRUE,
        observed = observed
      )
    })
    if (!length(group_data)) {
      group_data <- list(list(
        values = rep(NA_real_, length(unit_ids)),
        complete = FALSE,
        eligible = FALSE,
        observed = NA_real_
      ))
    }
    group_matrix <- do.call(rbind, lapply(group_data, `[[`, "values"))
    complete <- vapply(group_data, `[[`, logical(1), "complete")
    eligible <- vapply(group_data, `[[`, logical(1), "eligible")
    observed <- vapply(group_data, `[[`, numeric(1), "observed")
    values <- rep(NA_real_, length(unit_ids))
    if (any(eligible)) {
      selected <- group_matrix[eligible, , drop = FALSE]
      values <- if (identical(or_method, "max")) {
        apply(selected, 2L, function(x) {
          x <- x[is.finite(x)]
          if (!length(x)) NA_real_ else max(x)
        })
      } else if (identical(or_method, "prob_or")) {
        apply(selected, 2L, function(x) {
          x <- x[is.finite(x)]
          if (!length(x)) return(NA_real_)
          1 - prod(1 - pmin(pmax(x, 0), 1))
        })
      } else {
        apply(selected, 2L, function(x) {
          x <- x[is.finite(x)]
          if (!length(x)) return(NA_real_)
          sum(x) / sqrt(length(x))
        })
      }
    }
    all_genes <- unique(tolower(unlist(groups, use.names = FALSE)))
    all_genes <- all_genes[!is.na(all_genes) & nzchar(all_genes)]
    observed_genes <- intersect(all_genes, rownames(gene_confidence))
    coverage <- if (!length(all_genes)) {
      NA_real_
    } else {
      length(observed_genes) / length(all_genes)
    }
    data.frame(
      reaction_id = rid,
      pool_id = unit_ids,
      reaction_confidence = values,
      reaction_evidence_score = values,
      reaction_confidence_unpenalized = values,
      reaction_evidence_score_unpenalized = values,
      confidence_source = "gpr_aware_gene_confidence",
      evidence_source = "gpr_aware_gene_confidence",
      n_gpr_genes_total = length(all_genes),
      n_gpr_genes_multiome = length(observed_genes),
      multiome_coverage_fraction = coverage,
      missing_gpr_gene_fraction = if (is.na(coverage)) NA_real_ else 1 - coverage,
      missing_gene_fraction = if (is.na(coverage)) NA_real_ else 1 - coverage,
      missing_subunit_confidence_penalty = NA_real_,
      detection_available = !is.null(unit_detection),
      mean_gpr_detection_rate = NA_real_,
      low_confidence_threshold = if (is.null(low_confidence_threshold)) {
        NA_real_
      } else {
        low_confidence_threshold
      },
      n_and_groups_total = length(groups),
      n_and_groups_complete = sum(complete),
      n_and_groups_eligible = sum(eligible),
      complete_and_group_fraction = if (!length(groups)) NA_real_ else mean(complete),
      best_and_group_observed_fraction = if (all(is.na(observed))) {
        NA_real_
      } else {
        max(observed, na.rm = TRUE)
      },
      any_incomplete_gpr_group_flag = any(!complete),
      reaction_unsupported_by_complete_gpr_flag = !any(complete),
      missing_required_subunit_flag = any(!complete),
      no_complete_gpr_group_flag = !any(complete),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (!is.null(low_confidence_threshold)) {
    out$low_confidence_reaction_flag <-
      is.finite(out$reaction_confidence) &
      out$reaction_confidence < low_confidence_threshold
    return(out)
  }
  rc_apply_low_confidence_quantile(out, low_confidence_quantile)
}

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
      low_confidence_quantile = low_confidence_quantile,
      low_confidence_threshold = low_confidence_threshold
    )
    return(rc_set_confidence_source(out, source, detection_available = detection_available))
  }
}

rc_gpr_gene_ids <- function(gpr_list) {
  genes <- unique(tolower(unlist(gpr_list, use.names = FALSE)))
  genes[!is.na(genes) & nzchar(genes)]
}

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

