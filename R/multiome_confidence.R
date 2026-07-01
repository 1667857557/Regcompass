#' Stratum-wise empirical percentiles for feature-by-pool matrices
#'
#' Percentiles are computed within each requested pool stratum so cell-type or
#' condition composition does not dominate RNA/ATAC confidence components.
#' @export
rc_percentile_by_stratum <- function(x, pool_meta = NULL, stratum_col = NULL) {
  x <- as.matrix(x)
  if (is.null(colnames(x))) stop("`x` must have pool IDs in colnames().", call. = FALSE)
  if (is.null(stratum_col)) {
    strata <- rep("global", ncol(x))
  } else {
    if (is.null(pool_meta) || !"pool_id" %in% colnames(pool_meta) || !stratum_col %in% colnames(pool_meta)) {
      stop("`pool_meta` must contain `pool_id` and `stratum_col`.", call. = FALSE)
    }
    pool_meta <- pool_meta[match(colnames(x), pool_meta$pool_id), , drop = FALSE]
    if (anyNA(pool_meta$pool_id)) stop("`pool_meta` is missing metadata for some matrix columns.", call. = FALSE)
    strata <- as.character(pool_meta[[stratum_col]])
  }
  out <- x
  for (st in unique(strata)) {
    cols <- which(strata == st)
    out[, cols] <- t(apply(x[, cols, drop = FALSE], 1, rc_percentile_vector))
  }
  out
}

rc_percentile_vector <- function(v) {
  ok <- is.finite(v)
  out <- rep(NA_real_, length(v))
  n <- sum(ok)
  if (n == 0L) return(out)
  if (n == 1L) {
    out[ok] <- 1
    return(out)
  }
  out[ok] <- (rank(v[ok], ties.method = "average") - 1) / (n - 1)
  out
}

#' RNA-ATAC concordance with discrete-rank null correction
#' @export
rc_concordance_null_correct <- function(p_rna, p_atac) {
  p_rna <- as.matrix(p_rna); p_atac <- as.matrix(p_atac)
  if (!identical(dim(p_rna), dim(p_atac))) stop("`p_rna` and `p_atac` must have identical dimensions.", call. = FALSE)
  concord <- 1 - abs(p_rna - p_atac)
  n <- rowSums(is.finite(p_rna) & is.finite(p_atac))
  e_null <- 2 / 3 + 1 / (3 * n^2)
  e_null[!is.finite(e_null) | n < 1L] <- NA_real_
  denom <- 1 - e_null
  out <- sweep(sweep(concord, 1, e_null, "-"), 1, denom, "/")
  out <- pmax(0, pmin(1, out))
  out[!is.finite(out)] <- NA_real_
  out
}

#' Fisher-z shrinkage for row-wise Spearman correlations
#' @export
rc_fisher_shrink <- function(x, y, n0 = 30) {
  x <- as.matrix(x); y <- as.matrix(y)
  if (!identical(dim(x), dim(y))) stop("`x` and `y` must have identical dimensions.", call. = FALSE)
  res <- t(vapply(seq_len(nrow(x)), function(i) {
    ok <- is.finite(x[i, ]) & is.finite(y[i, ])
    n <- sum(ok)
    if (n < 4L) return(c(rho = NA_real_, rho_shrink = 0, n = n, low_correlation_power_flag = n < 10L))
    rho <- suppressWarnings(stats::cor(x[i, ok], y[i, ok], method = "spearman"))
    if (!is.finite(rho)) rho <- 0
    rho_clip <- min(0.999, max(-0.999, rho))
    lambda <- max(0, (n - 3) / (n - 3 + n0))
    rho_shrink <- tanh(lambda * atanh(rho_clip))
    c(rho = rho, rho_shrink = rho_shrink, n = n, low_correlation_power_flag = n < 10L)
  }, numeric(4)))
  out <- as.data.frame(res)
  out$low_correlation_power_flag <- as.logical(out$low_correlation_power_flag)
  rownames(out) <- rownames(x)
  out
}

#' Positive ATAC link confidence from normalized peak accessibility percentiles
#' @export
rc_link_confidence <- function(p_atac_peak, peak_gene_links) {
  p_atac_peak <- as.matrix(p_atac_peak)
  required <- c("peak_id", "gene", "weight")
  missing <- setdiff(required, colnames(peak_gene_links))
  if (length(missing) > 0L) stop("`peak_gene_links` is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  genes <- unique(as.character(peak_gene_links$gene))
  out <- matrix(0, nrow = length(genes), ncol = ncol(p_atac_peak), dimnames = list(genes, colnames(p_atac_peak)))
  for (g in genes) {
    links <- peak_gene_links[as.character(peak_gene_links$gene) == g, , drop = FALSE]
    links <- links[links$peak_id %in% rownames(p_atac_peak), , drop = FALSE]
    w <- pmax(0, as.numeric(links$weight))
    if (nrow(links) == 0L || sum(w) <= 0) next
    w <- w / sum(w)
    out[g, ] <- as.numeric(crossprod(w, p_atac_peak[links$peak_id, , drop = FALSE]))
  }
  attr(out, "repressive_link_flag") <- any(as.numeric(peak_gene_links$weight) < 0, na.rm = TRUE)
  out
}

#' Nonnegative gene-level multiome confidence
#' @export
rc_gene_confidence <- function(concord_ra_norm,
                               rel_ra_pos,
                               det_rna,
                               link_conf = NULL,
                               qc = NULL,
                               gpr_gene_observed = NULL,
                               concord_rt_norm = NULL,
                               rel_rt_pos = NULL) {
  base <- as.matrix(concord_ra_norm)
  det_rna <- as.matrix(det_rna)
  if (!identical(dim(base), dim(det_rna))) stop("`concord_ra_norm` and `det_rna` must have identical dimensions.", call. = FALSE)
  zeros <- matrix(0, nrow = nrow(base), ncol = ncol(base), dimnames = dimnames(base))
  if (is.null(link_conf)) link_conf <- zeros else link_conf <- as.matrix(link_conf)
  if (is.null(concord_rt_norm)) concord_rt_norm <- zeros else concord_rt_norm <- as.matrix(concord_rt_norm)
  if (is.null(rel_rt_pos)) rel_rt_pos <- rep(0, nrow(base))
  if (is.null(gpr_gene_observed)) gpr_gene_observed <- matrix(1, nrow = nrow(base), ncol = ncol(base), dimnames = dimnames(base))
  if (is.null(qc)) qc <- rep(1, ncol(base))
  if (!identical(dim(base), dim(link_conf)) || !identical(dim(base), dim(concord_rt_norm)) || !identical(dim(base), dim(as.matrix(gpr_gene_observed)))) {
    stop("All matrix confidence components must have identical dimensions.", call. = FALSE)
  }
  rel_ra_pos <- pmax(0, as.numeric(rel_ra_pos))
  rel_rt_pos <- pmax(0, as.numeric(rel_rt_pos))
  if (length(rel_ra_pos) != nrow(base) || length(rel_rt_pos) != nrow(base)) stop("Reliability vectors must have one value per gene/row.", call. = FALSE)
  if (length(qc) != ncol(base)) stop("`qc` must have one value per pool/column.", call. = FALSE)
  qc_mat <- matrix(pmax(0, pmin(1, as.numeric(qc))), nrow = nrow(base), ncol = ncol(base), byrow = TRUE)
  conf <- 0.25 * sweep(pmax(0, pmin(1, base)), 1, rel_ra_pos, "*") +
    0.15 * sweep(pmax(0, pmin(1, concord_rt_norm)), 1, rel_rt_pos, "*") +
    0.20 * pmax(0, pmin(1, det_rna)) +
    0.15 * pmax(0, pmin(1, link_conf)) +
    0.15 * qc_mat +
    0.10 * pmax(0, pmin(1, as.matrix(gpr_gene_observed)))
  pmax(0, pmin(1, conf))
}
