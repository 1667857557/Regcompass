#' Continuous shrinkage Q95 calibration
#' @export
rc_q95_shrink <- function(C_raw, pool_meta = NULL, stratum_col = NULL, q = 0.95, n0 = 80, eps = 1e-6) {
  C_raw <- as.matrix(C_raw)
  if (is.null(rownames(C_raw)) || is.null(colnames(C_raw))) stop("`C_raw` must have reaction rownames and pool colnames.", call. = FALSE)
  global_q <- apply(C_raw, 1, stats::quantile, probs = q, na.rm = TRUE, names = FALSE)
  if (!is.null(stratum_col)) {
    if (is.null(pool_meta) || !stratum_col %in% colnames(pool_meta)) stop("`pool_meta` with `stratum_col` is required.", call. = FALSE)
    pool_meta <- pool_meta[match(colnames(C_raw), pool_meta$pool_id), , drop = FALSE]
    strata <- as.character(pool_meta[[stratum_col]])
  } else {
    strata <- rep("global", ncol(C_raw))
  }
  diag_list <- list(); C_rel <- C_raw
  idx <- 1L
  for (st in unique(strata)) {
    pools <- which(strata == st)
    n <- length(pools); rho <- n / (n + n0)
    q_st <- apply(C_raw[, pools, drop = FALSE], 1, stats::quantile, probs = q, na.rm = TRUE, names = FALSE)
    q_shrink <- rho * q_st + (1 - rho) * global_q
    C_rel[, pools] <- sweep(C_raw[, pools, drop = FALSE], 1, q_shrink + eps, "/")
    diag_list[[idx]] <- data.frame(reaction_id = rownames(C_raw), stratum = st, n = n,
      q_stratum = as.numeric(q_st), q_global = as.numeric(global_q), rho_n = rho,
      q_shrink = as.numeric(q_shrink), q95_very_low_power = n < 5L,
      q95_low_power = n < 20L, q95_moderate_power = n < 100L,
      q95_high_power = n >= 400L, stringsAsFactors = FALSE)
    idx <- idx + 1L
  }
  C_rel[C_rel > 1] <- 1
  list(C_rel = C_rel, Q = do.call(rbind, diag_list))
}

#' Calibrate raw reaction capacities by continuous reaction-wise Q95 shrinkage
#' @export
rc_q95_calibrate <- function(C_raw, min_direct = 100, eps = 1e-6, bootstrap = TRUE, B = 200, BPPARAM = NULL, n0 = 80) {
  C_raw <- as.matrix(C_raw)
  out <- rc_q95_shrink(C_raw, q = 0.95, n0 = n0, eps = eps)
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
rc_q95_bootstrap <- function(x, B = 200) {
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
                                   promiscuity_mode = c("sqrt", "linear", "none"),
                                   tau = 0.08,
                                   min_direct = 100,
                                   bootstrap = TRUE,
                                   B = 200,
                                   BPPARAM = NULL) {
  promiscuity_mode <- match.arg(promiscuity_mode)
  if (is.null(rownames(pool_expression))) {
    stop("`pool_expression` must have gene IDs in rownames().", call. = FALSE)
  }
  parsed <- rc_parse_gpr_table(gpr_table)
  gene_score <- rc_sigmoid(rc_robust_z(pool_expression))
  rownames(gene_score) <- tolower(rownames(gene_score))

  C_raw <- rc_reaction_capacity(parsed, gene_score, promiscuity_mode = promiscuity_mode, tau = tau, BPPARAM = BPPARAM)
  calibrated <- rc_q95_calibrate(C_raw, min_direct = min_direct, bootstrap = bootstrap, B = B, BPPARAM = BPPARAM)
  gpr_diag <- rc_gpr_diagnostics(parsed, rownames(gene_score))
  confidence <- rc_reaction_confidence(parsed, pool_detection)

  list(
    reaction_capacity_L1 = calibrated$C_rel,
    reaction_capacity_raw = C_raw,
    reaction_confidence = confidence,
    q95_diagnostics = calibrated$Q,
    gpr_diagnostics = gpr_diag,
    parsed_gpr = parsed
  )
}

#' Compute reaction confidence from GPR gene detection rates
#'
#' @param gpr_list A named list of parsed GPR rules.
#' @param pool_detection Optional gene-by-pool detection-rate matrix.
#'
#' @return A reaction-level data.frame. If `pool_detection` is supplied, includes
#' mean detection across available GPR genes and pools.
#' @export
rc_reaction_confidence <- function(gpr_list, pool_detection = NULL) {
  if (is.null(pool_detection)) {
    all_genes <- unique(unlist(gpr_list, use.names = FALSE))
    diag <- rc_gpr_diagnostics(gpr_list, all_genes)
    diag$detection_available <- FALSE
    diag$mean_gpr_detection_rate <- NA_real_
    return(diag)
  }
  pool_detection <- as.matrix(pool_detection)
  if (is.null(rownames(pool_detection))) {
    stop("`pool_detection` must have gene IDs in rownames().", call. = FALSE)
  }
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
