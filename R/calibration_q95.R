#' Continuous shrinkage Q95 calibration
#' @export
rc_q95_shrink <- function(C_raw, pool_meta = NULL, stratum_col = NULL, q = 0.95, n0 = 80, eps = 1e-6) {
  C_raw <- as.matrix(C_raw)
  if (is.null(rownames(C_raw)) || is.null(colnames(C_raw))) stop("`C_raw` must have reaction rownames and pool colnames.", call. = FALSE)
  global_q <- apply(C_raw, 1, rc_safe_quantile, probs = q)
  global_n <- rowSums(is.finite(C_raw))
  if (!is.null(stratum_col)) {
    if (is.null(pool_meta) || !"pool_id" %in% colnames(pool_meta) || !stratum_col %in% colnames(pool_meta)) stop("`pool_meta` with `pool_id` and `stratum_col` is required.", call. = FALSE)
    pool_meta <- pool_meta[match(colnames(C_raw), pool_meta$pool_id), , drop = FALSE]
    if (anyNA(pool_meta$pool_id)) stop("`pool_meta` is missing metadata for some capacity columns.", call. = FALSE)
    strata <- as.character(pool_meta[[stratum_col]])
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
      q95_very_low_power = n < 5L, q95_low_power = n < 20L,
      q95_moderate_power = n < 100L, q95_high_power = n >= 400L, stringsAsFactors = FALSE)
    idx <- idx + 1L
  }
  C_rel[C_rel > 1] <- 1
  list(C_rel = C_rel, Q = do.call(rbind, diag_list))
}

#' Calibrate raw reaction capacities by continuous reaction-wise Q95 shrinkage
#' @export
rc_q95_calibrate <- function(C_raw, min_direct = 100, eps = 1e-6, bootstrap = TRUE, B = 500, BPPARAM = NULL, n0 = 80, pool_meta = NULL, stratum_col = NULL) {
  C_raw <- as.matrix(C_raw)
  out <- rc_q95_shrink(C_raw, pool_meta = pool_meta, stratum_col = stratum_col, q = 0.95, n0 = n0, eps = eps)
  Q <- out$Q
  names(Q)[names(Q) == "q_shrink"] <- "q_value"
  Q$quantile_used <- 0.95
  Q$n_finite <- rowSums(is.finite(C_raw))[match(Q$reaction_id, rownames(C_raw))]
  Q$low_n_flag <- Q$n_finite < 20L
  if (isTRUE(bootstrap)) {
    boot <- rc_q95_bootstrap_diagnostics(C_raw, Q, pool_meta = pool_meta, stratum_col = stratum_col, B = B, BPPARAM = BPPARAM)
    Q$q95_bootstrap <- boot[, "q95"]; Q$q95_ci_low <- boot[, "ci_low"]
    Q$q95_ci_high <- boot[, "ci_high"]; Q$q95_ci_width <- boot[, "width"]
    Q$q95_unstable_flag <- is.finite(Q$q95_ci_width) & Q$q95_ci_width > Q$q_value
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

#' Compute reaction confidence from gene confidence or RNA detection
#'
#' Reaction confidence is the median available GPR-gene confidence per pool,
#' multiplied by the observed fraction of GPR genes. This preserves the capacity
#' estimate from available genes but downgrades reactions with missing GPR genes.
#' @export
rc_run_layer1_capacity <- function(gpr_table,
                                   pool_expression,
                                   pool_detection = NULL,
                                   pool_meta = NULL,
                                   stratum_col = NULL,
                                   gene_confidence = NULL,
                                   gene_discordance = NULL,
                                   and_methods = c("min", "boltzmann_0.08", "boltzmann_0.20", "mean"),
                                   tau_sensitivity = c(0.08, 0.20),
                                   promiscuity_sensitivity = c("none", "sqrt", "linear"),
                                   promiscuity_mode = c("sqrt", "linear", "none"),
                                   tau = 0.20,
                                   and_method = c("boltzmann", "min", "mean"),
                                   min_direct = 100,
                                   bootstrap = TRUE,
                                   B = 500,
                                   BPPARAM = NULL) {
  promiscuity_mode <- match.arg(promiscuity_mode)
  and_method <- match.arg(and_method)
  if (is.null(rownames(pool_expression))) {
    stop("`pool_expression` must have gene IDs in rownames().", call. = FALSE)
  }
  gene_confidence <- as.matrix(gene_confidence)
  if (is.null(rownames(gene_confidence)) || is.null(colnames(gene_confidence))) stop("Gene confidence/detection must have gene rownames and pool colnames.", call. = FALSE)
  rownames(gene_confidence) <- tolower(rownames(gene_confidence))

  rows <- lapply(names(gpr_list), function(rid) {
    genes_all <- unique(unlist(gpr_list[[rid]], use.names = FALSE))
    genes_all <- genes_all[nzchar(genes_all)]
    total <- length(genes_all)
    genes <- intersect(genes_all, rownames(gene_confidence))
    missing_frac <- if (total == 0L) NA_real_ else 1 - length(genes) / total
    raw <- if (length(genes) == 0L) rep(NA_real_, ncol(gene_confidence)) else matrixStats::colMedians(gene_confidence[genes, , drop = FALSE], na.rm = TRUE)
    conf <- raw * ifelse(is.na(missing_frac), NA_real_, 1 - missing_frac)
    data.frame(
      reaction_id = rid,
      pool_id = colnames(gene_confidence),
      reaction_confidence = pmax(0, pmin(1, conf)),
      raw_reaction_confidence = pmax(0, pmin(1, raw)),
      missing_gpr_gene_fraction = missing_frac,
      low_confidence_reaction_flag = conf < 0.25 | is.na(conf),
      confidence_source = source,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Run the simplified RegCompassR Layer 1 workflow from RNA raw counts
#' @export
rc_run_layer1_from_counts <- function(gpr_table,
                                      rna_counts,
                                      pool_map,
                                      pool_meta = NULL,
                                      stratum_col = "cell_type",
                                      tau = 0.20,
                                      n0 = 80,
                                      BPPARAM = NULL) {
  pb <- rc_pseudobulk_counts(rna_counts, pool_map, fun = "sum", BPPARAM = BPPARAM)
  if (is.null(pool_meta)) pool_meta <- rc_build_pool_metadata(pool_map)
  filtered <- rc_filter_empty_pools(pb, pool_meta)
  rna_logcpm <- rc_logcpm(filtered$counts)
  pool_meta <- filtered$pool_meta

  detection <- rc_pool_detection(rna_counts, pool_map, BPPARAM = BPPARAM)
  detection <- detection[, colnames(rna_logcpm), drop = FALSE]

  parsed <- rc_parse_gpr_table(gpr_table)
  gene_score <- rc_gene_score(pool_expression)
  rownames(gene_score) <- tolower(rownames(gene_score))

  C_raw <- rc_reaction_capacity(parsed, gene_score, promiscuity_mode = promiscuity_mode, tau = tau, and_method = and_method, BPPARAM = BPPARAM)
  calibrated <- rc_q95_calibrate(C_raw, min_direct = min_direct, bootstrap = bootstrap, B = B, BPPARAM = BPPARAM, pool_meta = pool_meta, stratum_col = stratum_col)
  gpr_diag <- rc_gpr_diagnostics(parsed, rownames(gene_score))
  confidence <- rc_reaction_confidence(parsed, gene_confidence = gene_confidence, pool_detection = pool_detection)
  tau_sens <- rc_capacity_sensitivity(parsed, gene_score, variable = "tau", values = tau_sensitivity, promiscuity_mode = promiscuity_mode, and_method = and_method, BPPARAM = BPPARAM)
  prom_sens <- rc_capacity_sensitivity(parsed, gene_score, variable = "promiscuity", values = promiscuity_sensitivity, tau = tau, and_method = and_method, BPPARAM = BPPARAM)
  and_long <- rc_and_method_capacity_long(parsed, gene_score, and_methods = and_methods, promiscuity_mode = promiscuity_mode, BPPARAM = BPPARAM)
  and_sens <- rc_and_method_sensitivity(and_long)

  list(
    C_raw = C_raw,
    C_rel = calibrated$C_rel,
    reaction_confidence = confidence,
    capacity_long = and_long,
    and_method_sensitivity = and_sens,
    q95_diagnostics = calibrated$Q,
    gpr_diagnostics = gpr_diag,
    gene_discordance = gene_discordance,
    tau_sensitivity = tau_sens,
    promiscuity_sensitivity = prom_sens,
    parsed_gpr = parsed
  )
}

#' Compute reaction confidence from gene confidence or detection fallback
#' @export
rc_reaction_confidence <- function(gpr_list, gene_confidence = NULL, pool_detection = NULL) {
  if (!is.null(gene_confidence)) {
    gene_confidence <- as.matrix(gene_confidence)
    rownames(gene_confidence) <- tolower(rownames(gene_confidence))
    rows <- lapply(names(gpr_list), function(rid) {
      genes <- unique(unlist(gpr_list[[rid]], use.names = FALSE))
      total <- length(genes)
      genes <- intersect(genes, rownames(gene_confidence))
      vals <- if (length(genes) == 0L) rep(NA_real_, ncol(gene_confidence)) else matrixStats::colMedians(gene_confidence[genes, , drop = FALSE], na.rm = TRUE)
      miss <- if (total == 0L) NA_real_ else 1 - length(genes) / total
      penalty <- ifelse(is.na(miss), NA_real_, miss)
      vals_penalized <- vals * (1 - ifelse(is.na(penalty), 0, penalty))
      data.frame(reaction_id = rid, pool_id = colnames(gene_confidence), reaction_confidence = vals_penalized,
                 reaction_confidence_unpenalized = vals, missing_gpr_gene_fraction = miss,
                 missing_subunit_confidence_penalty = penalty,
                 low_confidence_reaction_flag = vals_penalized < 0.25, stringsAsFactors = FALSE)
    })
    return(do.call(rbind, rows))
  }
  if (is.null(pool_detection)) {
    all_genes <- unique(unlist(gpr_list, use.names = FALSE))
    diag <- rc_gpr_diagnostics(gpr_list, all_genes)
    diag$pool_id <- NA_character_
    diag$reaction_confidence <- NA_real_
    diag$detection_available <- FALSE
    diag$mean_gpr_detection_rate <- NA_real_
    return(diag)
  }
  pool_detection <- as.matrix(pool_detection)
  if (is.null(rownames(pool_detection))) stop("`pool_detection` must have gene IDs in rownames().", call. = FALSE)
  rownames(pool_detection) <- tolower(rownames(pool_detection))
  rows <- lapply(names(gpr_list), function(rid) {
    genes_all <- unique(unlist(gpr_list[[rid]], use.names = FALSE))
    total <- length(genes_all)
    genes <- intersect(genes_all, rownames(pool_detection))
    vals <- if (length(genes) == 0L) rep(NA_real_, ncol(pool_detection)) else matrixStats::colMedians(pool_detection[genes, , drop = FALSE], na.rm = TRUE)
    miss <- if (total == 0L) NA_real_ else 1 - length(genes) / total
    penalty <- ifelse(is.na(miss), NA_real_, miss)
    vals_penalized <- vals * (1 - ifelse(is.na(penalty), 0, penalty))
    data.frame(reaction_id = rid, pool_id = colnames(pool_detection), reaction_confidence = vals_penalized,
               reaction_confidence_unpenalized = vals, missing_gpr_gene_fraction = miss,
               missing_subunit_confidence_penalty = penalty,
               detection_available = TRUE, mean_gpr_detection_rate = vals, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Compute compact capacity sensitivity summaries
#' @export
rc_capacity_sensitivity <- function(gpr_list, gene_score, variable = c("tau", "promiscuity"), values, promiscuity_mode = "sqrt", tau = 0.20, and_method = "boltzmann", BPPARAM = NULL) {
  variable <- match.arg(variable)
  if (length(values) == 0L) return(data.frame())
  rows <- lapply(values, function(v) {
    C <- if (variable == "tau") {
      rc_reaction_capacity(gpr_list, gene_score, promiscuity_mode = promiscuity_mode, tau = as.numeric(v), and_method = and_method, BPPARAM = BPPARAM)
    } else {
      rc_reaction_capacity(gpr_list, gene_score, promiscuity_mode = as.character(v), tau = tau, and_method = and_method, BPPARAM = BPPARAM)
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


rc_q95_bootstrap_diagnostics <- function(C_raw, Q, pool_meta = NULL, stratum_col = NULL, B = 500, BPPARAM = NULL) {
  C_raw <- as.matrix(C_raw)
  strata <- if (!is.null(stratum_col)) {
    pool_meta <- pool_meta[match(colnames(C_raw), pool_meta$pool_id), , drop = FALSE]
    as.character(pool_meta[[stratum_col]])
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
rc_and_method_capacity_long <- function(gpr_list, gene_score, and_methods = c("min", "boltzmann_0.08", "boltzmann_0.20", "mean"), promiscuity_mode = "sqrt", BPPARAM = NULL) {
  rows <- lapply(and_methods, function(m) {
    spec <- rc_parse_and_method(m)
    C <- rc_reaction_capacity(gpr_list, gene_score, promiscuity_mode = promiscuity_mode, tau = spec$tau, and_method = spec$method, BPPARAM = BPPARAM)
    data.frame(reaction_id = rep(rownames(C), times = ncol(C)), pool_id = rep(colnames(C), each = nrow(C)), and_method = m, C_raw = as.vector(C), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Summarize AND-method sensitivity by reaction and pool
#' @export
rc_and_method_sensitivity <- function(capacity_long) {
  split_key <- interaction(capacity_long$reaction_id, capacity_long$pool_id, drop = TRUE)
  rows <- lapply(split(capacity_long, split_key), function(df) {
    tau_vals <- df$C_raw[df$and_method %in% c("boltzmann_0.08", "boltzmann_0.20")]
    data.frame(reaction_id = df$reaction_id[[1]], pool_id = df$pool_id[[1]],
               capacity_range = max(df$C_raw, na.rm = TRUE) - min(df$C_raw, na.rm = TRUE),
               tau_sensitive_flag = length(tau_vals) >= 2L && diff(range(tau_vals, na.rm = TRUE)) > 0.05,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

rc_parse_and_method <- function(x) {
  if (x == "min") return(list(method = "min", tau = 0.20))
  if (x == "mean") return(list(method = "mean", tau = 0.20))
  if (grepl("^boltzmann_", x)) return(list(method = "boltzmann", tau = as.numeric(sub("^boltzmann_", "", x))))
  stop("Unsupported AND method: ", x, call. = FALSE)
}

#' End-to-end Layer 1 run from raw counts and a pool map
#' @export
rc_run_layer1_from_counts <- function(gpr_table,
                                      rna_counts,
                                      pool_map,
                                      pool_meta = NULL,
                                      atac_counts = NULL,
                                      peak_gene_links = NULL,
                                      stratum_col = "cell_type",
                                      promiscuity_mode = "sqrt",
                                      and_method = "boltzmann",
                                      tau = 0.20,
                                      bootstrap = TRUE,
                                      B = 500,
                                      BPPARAM = NULL) {
  rna_pb <- rc_pseudobulk_counts(rna_counts, pool_map, fun = "sum", BPPARAM = BPPARAM)
  if (is.null(pool_meta)) pool_meta <- rc_build_pool_metadata(pool_map)
  filtered <- rc_filter_empty_pools(rna_pb, pool_meta)
  rna_logcpm <- rc_logcpm(filtered$counts)
  pool_meta <- filtered$pool_meta
  rna_detection <- rc_pool_detection(rna_counts, pool_map, BPPARAM = BPPARAM)
  rna_detection <- rna_detection[, colnames(rna_logcpm), drop = FALSE]

  gene_conf <- NULL
  confidence_source <- "rna_detection"
  if (!is.null(atac_counts) && !is.null(peak_gene_links)) {
    atac_peak <- rc_atac_pool_logcpm(atac_counts, pool_map, min_pools = 3, BPPARAM = BPPARAM)
    common_pools <- intersect(colnames(rna_logcpm), colnames(atac_peak))
    rna_logcpm <- rna_logcpm[, common_pools, drop = FALSE]
    rna_detection <- rna_detection[, common_pools, drop = FALSE]
    pool_meta <- pool_meta[match(common_pools, pool_meta$pool_id), , drop = FALSE]
    p_rna <- rc_percentile_by_stratum(rna_logcpm, pool_meta = pool_meta, stratum_col = stratum_col)
    p_atac_peak <- rc_percentile_by_stratum(atac_peak[, common_pools, drop = FALSE], pool_meta = pool_meta, stratum_col = stratum_col)
    link_conf <- rc_link_confidence(p_atac_peak, peak_gene_links)
    genes <- intersect(rownames(p_rna), rownames(link_conf))
    if (length(genes) > 0L) {
      concord <- rc_concordance_null_correct(p_rna[genes, , drop = FALSE], link_conf[genes, , drop = FALSE], pool_meta = pool_meta, stratum_col = stratum_col)
      rel <- rc_fisher_shrink(rna_logcpm[genes, , drop = FALSE], link_conf[genes, , drop = FALSE])$rel_positive
      names(rel) <- genes
      gene_conf <- rc_gene_confidence(concord, rel_ra_pos = rel, det_rna = rna_detection[genes, , drop = FALSE], link_conf = link_conf[genes, , drop = FALSE])
      confidence_source <- "multiome_link_confidence"
    }
  }

  out <- rc_run_layer1_capacity(gpr_table = gpr_table, pool_expression = rna_logcpm, pool_detection = rna_detection,
                                pool_meta = pool_meta, stratum_col = stratum_col, gene_confidence = gene_conf,
                                promiscuity_mode = promiscuity_mode, and_method = and_method, tau = tau, bootstrap = bootstrap, B = B, BPPARAM = BPPARAM)
  out$pool_meta <- pool_meta
  out$reaction_confidence_source <- confidence_source
  out
}
