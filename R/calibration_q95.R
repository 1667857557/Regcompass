#' Calibrate raw reaction capacities by reaction-wise upper quantiles
#'
#' Reactions with at least `min_direct` finite pool values use Q95; smaller
#' reactions use Q90 as a conservative fallback for early v0.3 diagnostics.
#'
#' @param C_raw Reaction-by-pool raw capacity matrix.
#' @param min_direct Minimum number of finite values needed for direct Q95.
#' @param eps Positive stabilizer added to the quantile denominator.
#'
#' @return A list with `C_rel`, the clipped relative capacity matrix, and `Q`, a
#' data.frame with reaction-wise quantiles and diagnostic metadata.
#' @export
rc_q95_calibrate <- function(C_raw, min_direct = 100, eps = 1e-6) {
  C_raw <- as.matrix(C_raw)
  if (is.null(rownames(C_raw))) {
    stop("`C_raw` must have reaction IDs in rownames().", call. = FALSE)
  }
  if (!is.numeric(min_direct) || length(min_direct) != 1L || is.na(min_direct) || min_direct < 1) {
    stop("`min_direct` must be a single positive number.", call. = FALSE)
  }

  q <- apply(C_raw, 1, function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
      return(NA_real_)
    }
    if (length(x) >= min_direct) {
      stats::quantile(x, 0.95, na.rm = TRUE, names = FALSE)
    } else {
      stats::quantile(x, 0.90, na.rm = TRUE, names = FALSE)
    }
  })

  n_finite <- rowSums(is.finite(C_raw))
  quantile_used <- ifelse(n_finite >= min_direct, 0.95, 0.90)
  quantile_used[n_finite == 0L] <- NA_real_
  C_rel <- sweep(C_raw, 1, q + eps, "/")
  C_rel[C_rel > 1] <- 1

  Q <- data.frame(
    reaction_id = rownames(C_raw),
    q_value = as.numeric(q),
    quantile_used = quantile_used,
    n_finite = as.integer(n_finite),
    low_n_flag = n_finite < min_direct,
    stringsAsFactors = FALSE
  )
  list(C_rel = C_rel, Q = Q)
}

#' Run v0.3 Layer 1 GPR capacity workflow
#'
#' @param gpr_table A data.frame with `reaction_id` and `gpr` columns.
#' @param pool_expression Gene-by-pool expression score input. Use normalized
#' data or residuals, not imputed matrices.
#' @param pool_detection Optional gene-by-pool detection-rate matrix from raw
#' counts for confidence summaries.
#' @param promiscuity_mode One of `"sqrt"`, `"linear"`, or `"none"`.
#' @param tau Positive Boltzmann temperature for AND complexes.
#' @param min_direct Minimum finite pool count for Q95 calibration.
#' @param BPPARAM Optional `BiocParallelParam` used to parallelize reactions.
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
                                   BPPARAM = NULL) {
  promiscuity_mode <- match.arg(promiscuity_mode)
  if (is.null(rownames(pool_expression))) {
    stop("`pool_expression` must have gene IDs in rownames().", call. = FALSE)
  }
  parsed <- rc_parse_gpr_table(gpr_table)
  gene_score <- rc_sigmoid(rc_robust_z(pool_expression))
  rownames(gene_score) <- tolower(rownames(gene_score))

  C_raw <- rc_reaction_capacity(parsed, gene_score, promiscuity_mode = promiscuity_mode, tau = tau, BPPARAM = BPPARAM)
  calibrated <- rc_q95_calibrate(C_raw, min_direct = min_direct)
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
    return(rc_gpr_diagnostics(gpr_list, character(0)))
  }
  pool_detection <- as.matrix(pool_detection)
  if (is.null(rownames(pool_detection))) {
    stop("`pool_detection` must have gene IDs in rownames().", call. = FALSE)
  }
  rownames(pool_detection) <- tolower(rownames(pool_detection))
  diag <- rc_gpr_diagnostics(gpr_list, rownames(pool_detection))
  mean_detection <- vapply(names(gpr_list), function(rid) {
    genes <- unique(unlist(gpr_list[[rid]], use.names = FALSE))
    genes <- intersect(genes, rownames(pool_detection))
    if (length(genes) == 0L) return(NA_real_)
    mean(pool_detection[genes, , drop = FALSE], na.rm = TRUE)
  }, numeric(1))
  diag$mean_gpr_detection_rate <- mean_detection
  diag
}
