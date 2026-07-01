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
    q_shrink <- rho * q_st + (1 - rho) * global_q
    C_rel[, pools] <- sweep(C_raw[, pools, drop = FALSE], 1, q_shrink + eps, "/")
    diag_list[[idx]] <- data.frame(reaction_id = rownames(C_raw), stratum = st,
      n = as.integer(n), n_global = as.integer(global_n), q_stratum = as.numeric(q_st),
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
    boot <- rc_parallel_lapply(seq_len(nrow(C_raw)), function(i) rc_q95_bootstrap(C_raw[i, ], B = B), BPPARAM = BPPARAM)
    boot <- do.call(rbind, boot)
    Q$q95_bootstrap <- boot[, "q95"]; Q$q95_ci_low <- boot[, "ci_low"]
    Q$q95_ci_high <- boot[, "ci_high"]; Q$q95_ci_width <- boot[, "width"]
    Q$q95_unstable_flag <- is.finite(Q$q95_ci_width) & Q$q95_ci_width > Q$q_value
  }
  list(C_rel = out$C_rel, Q = Q)
}

#' Bootstrap Q95 confidence interval for one reaction
#'
#' @param x Numeric raw capacity values for one reaction across pools.
#' @param B Number of bootstrap resamples.
#'
#' @return Named numeric vector: `q95`, `ci_low`, `ci_high`, and `width`.
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

#' Run v0.4 Layer 1 GPR capacity workflow with diagnostics
#'
#' @param gpr_table A data.frame with `reaction_id` and `gpr` columns.
#' @param pool_expression Gene-by-pool expression score input. Use normalized
#' data or residuals, not imputed matrices.
#' @param pool_detection Optional gene-by-pool detection-rate matrix from raw
#' counts for confidence summaries.
#' @param promiscuity_mode One of `"sqrt"`, `"linear"`, or `"none"`.
#' @param tau Positive Boltzmann temperature for AND complexes.
#' @param min_direct Minimum finite pool count for Q95 calibration.
#' @param bootstrap Logical; if `TRUE`, include bootstrap Q95 confidence intervals.
#' @param B Number of bootstrap resamples when `bootstrap = TRUE`.
#' @param BPPARAM Optional `BiocParallelParam` used to parallelize reactions and bootstrap.
#'
#' @return A list with `reaction_capacity_L1`, `reaction_capacity_raw`,
#' `reaction_confidence`, `q95_diagnostics`, `gpr_diagnostics`, and `parsed_gpr`.
#' @export
rc_run_layer1_capacity <- function(gpr_table,
                                   pool_expression,
                                   pool_detection = NULL,
                                   pool_meta = NULL,
                                   stratum_col = NULL,
                                   gene_confidence = NULL,
                                   gene_discordance = NULL,
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
  parsed <- rc_parse_gpr_table(gpr_table)
  gene_score <- rc_gene_score(pool_expression)
  rownames(gene_score) <- tolower(rownames(gene_score))

  C_raw <- rc_reaction_capacity(parsed, gene_score, promiscuity_mode = promiscuity_mode, tau = tau, and_method = and_method, BPPARAM = BPPARAM)
  calibrated <- rc_q95_calibrate(C_raw, min_direct = min_direct, bootstrap = bootstrap, B = B, BPPARAM = BPPARAM, pool_meta = pool_meta, stratum_col = stratum_col)
  gpr_diag <- rc_gpr_diagnostics(parsed, rownames(gene_score))
  confidence <- rc_reaction_confidence(parsed, gene_confidence = gene_confidence, pool_detection = pool_detection)
  tau_sens <- rc_capacity_sensitivity(parsed, gene_score, variable = "tau", values = tau_sensitivity, promiscuity_mode = promiscuity_mode, and_method = and_method, BPPARAM = BPPARAM)
  prom_sens <- rc_capacity_sensitivity(parsed, gene_score, variable = "promiscuity", values = promiscuity_sensitivity, tau = tau, and_method = and_method, BPPARAM = BPPARAM)

  list(
    reaction_capacity_L1 = calibrated$C_rel,
    reaction_capacity_raw = C_raw,
    reaction_confidence = confidence,
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
      data.frame(reaction_id = rid, pool_id = colnames(gene_confidence), reaction_confidence = vals,
                 missing_gpr_gene_fraction = if (total == 0L) NA_real_ else 1 - length(genes) / total,
                 low_confidence_reaction_flag = vals < 0.25, stringsAsFactors = FALSE)
    })
    return(do.call(rbind, rows))
  }
  if (is.null(pool_detection)) {
    all_genes <- unique(unlist(gpr_list, use.names = FALSE))
    diag <- rc_gpr_diagnostics(gpr_list, all_genes)
    diag$detection_available <- FALSE
    diag$mean_gpr_detection_rate <- NA_real_
    return(diag)
  }
  pool_detection <- as.matrix(pool_detection)
  if (is.null(rownames(pool_detection))) stop("`pool_detection` must have gene IDs in rownames().", call. = FALSE)
  rownames(pool_detection) <- tolower(rownames(pool_detection))
  diag <- rc_gpr_diagnostics(gpr_list, rownames(pool_detection))
  diag$detection_available <- TRUE
  mean_detection <- vapply(names(gpr_list), function(rid) {
    genes <- unique(unlist(gpr_list[[rid]], use.names = FALSE))
    genes <- intersect(genes, rownames(pool_detection))
    if (length(genes) == 0L) return(NA_real_)
    mean(pool_detection[genes, , drop = FALSE], na.rm = TRUE)
  }, numeric(1))
  diag$mean_gpr_detection_rate <- mean_detection
  diag
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
