#' Continuous reaction-wise Q95 shrinkage calibration
#' @export
rc_q95_shrink <- function(C_raw, pool_meta = NULL, stratum_col = "cell_type", q = 0.95, n0 = 80, eps = 1e-6) {
  C_raw <- as.matrix(C_raw)
  if (is.null(rownames(C_raw)) || is.null(colnames(C_raw))) stop("`C_raw` must have reaction rownames and pool colnames.", call. = FALSE)
  global_q <- apply(C_raw, 1, rc_safe_quantile, probs = q)
  global_n <- rowSums(is.finite(C_raw))

  if (is.null(stratum_col)) {
    strata <- rep("global", ncol(C_raw))
  } else {
    if (is.null(pool_meta) || !"pool_id" %in% colnames(pool_meta) || !stratum_col %in% colnames(pool_meta)) {
      stop("`pool_meta` must contain `pool_id` and `", stratum_col, "` for cell-type Q95 calibration.", call. = FALSE)
    }
    pool_meta <- pool_meta[match(colnames(C_raw), pool_meta$pool_id), , drop = FALSE]
    if (anyNA(pool_meta$pool_id)) stop("`pool_meta` is missing metadata for some capacity columns.", call. = FALSE)
    strata <- as.character(pool_meta[[stratum_col]])
  }

  C_rel <- C_raw
  diag_list <- list()
  idx <- 1L
  for (st in unique(strata)) {
    pools <- which(strata == st)
    n <- rowSums(is.finite(C_raw[, pools, drop = FALSE]))
    rho <- n / (n + n0)
    q_st <- apply(C_raw[, pools, drop = FALSE], 1, rc_safe_quantile, probs = q)
    q_base <- ifelse(is.finite(q_st), q_st, global_q)
    q_shrink <- rho * q_base + (1 - rho) * global_q
    C_rel[, pools] <- sweep(C_raw[, pools, drop = FALSE], 1, q_shrink + eps, "/")
    diag_list[[idx]] <- data.frame(
      reaction_id = rownames(C_raw),
      stratum = st,
      n = as.integer(n),
      n_global = as.integer(global_n),
      q_stratum = as.numeric(q_st),
      q_global = as.numeric(global_q),
      rho_n = as.numeric(rho),
      q_shrink = as.numeric(q_shrink),
      q95_very_low_power = n < 5L,
      q95_low_power = n < 20L,
      q95_moderate_power = n < 100L,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }

  C_rel[C_rel > 1] <- 1
  list(C_rel = C_rel, Q = do.call(rbind, diag_list))
}

#' Calibrate raw reaction capacity by Q95 shrinkage
#' @export
rc_q95_calibrate <- function(C_raw, pool_meta = NULL, stratum_col = "cell_type", n0 = 80, eps = 1e-6) {
  out <- rc_q95_shrink(C_raw, pool_meta = pool_meta, stratum_col = stratum_col, q = 0.95, n0 = n0, eps = eps)
  names(out$Q)[names(out$Q) == "q_shrink"] <- "q_value"
  out$Q$quantile_used <- 0.95
  out
}

#' Compute reaction confidence from gene confidence or RNA detection
#'
#' Reaction confidence is the median available GPR-gene confidence per pool,
#' multiplied by the observed fraction of GPR genes. This preserves the capacity
#' estimate from available genes but downgrades reactions with missing GPR genes.
#' @export
rc_reaction_confidence <- function(gpr_list, gene_confidence = NULL, pool_detection = NULL) {
  source <- "gene_confidence"
  if (is.null(gene_confidence)) {
    if (is.null(pool_detection)) stop("Provide `gene_confidence` or `pool_detection`.", call. = FALSE)
    gene_confidence <- pool_detection
    source <- "rna_detection"
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
  gene_score <- rc_gene_score(rna_logcpm)
  rownames(gene_score) <- tolower(rownames(gene_score))

  C_raw <- rc_reaction_capacity(parsed, gene_score, tau = tau, BPPARAM = BPPARAM)
  calibrated <- rc_q95_calibrate(C_raw, pool_meta = pool_meta, stratum_col = stratum_col, n0 = n0)
  confidence <- rc_reaction_confidence(parsed, pool_detection = detection)
  gpr_diag <- rc_gpr_diagnostics(parsed, rownames(gene_score))

  C_min <- rc_hard_min_capacity(parsed, gene_score, BPPARAM = BPPARAM)
  tau_delta <- abs(C_raw - C_min)
  tau_sensitive <- data.frame(
    reaction_id = rownames(C_raw),
    tau = tau,
    mean_abs_delta_vs_hard_min = rowMeans(tau_delta, na.rm = TRUE),
    tau_sensitive_flag = rowMeans(tau_delta, na.rm = TRUE) > 0.15,
    stringsAsFactors = FALSE
  )

  pool_diag <- rc_pool_diagnostics(pool_map, rna_counts = rna_counts, gpr_genes = unique(unlist(parsed, use.names = FALSE)))

  list(
    C_raw = C_raw,
    C_rel = calibrated$C_rel,
    reaction_confidence = confidence,
    minimal_diagnostics = list(
      pool_diagnostics = pool_diag,
      q95_diagnostics = calibrated$Q,
      gpr_diagnostics = gpr_diag,
      tau_sensitivity = tau_sensitive
    ),
    pool_expression_logcpm = rna_logcpm,
    pool_detection = detection,
    pool_meta = pool_meta,
    parsed_gpr = parsed
  )
}

#' Backward-compatible Layer 1 function for already normalized pool expression
#' @export
rc_run_layer1_capacity <- function(gpr_table,
                                   pool_expression,
                                   pool_detection = NULL,
                                   pool_meta = NULL,
                                   stratum_col = "cell_type",
                                   tau = 0.20,
                                   n0 = 80,
                                   BPPARAM = NULL,
                                   ...) {
  parsed <- rc_parse_gpr_table(gpr_table)
  gene_score <- rc_gene_score(pool_expression)
  rownames(gene_score) <- tolower(rownames(gene_score))
  C_raw <- rc_reaction_capacity(parsed, gene_score, tau = tau, BPPARAM = BPPARAM)
  calibrated <- rc_q95_calibrate(C_raw, pool_meta = pool_meta, stratum_col = stratum_col, n0 = n0)
  confidence <- if (is.null(pool_detection)) NULL else rc_reaction_confidence(parsed, pool_detection = pool_detection)
  list(C_raw = C_raw, C_rel = calibrated$C_rel, reaction_confidence = confidence,
       minimal_diagnostics = list(q95_diagnostics = calibrated$Q, gpr_diagnostics = rc_gpr_diagnostics(parsed, rownames(gene_score))),
       parsed_gpr = parsed)
}

#' Layer 1 alias
#' @export
rc_layer1_capacity <- rc_run_layer1_from_counts

rc_safe_quantile <- function(x, probs) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
}
