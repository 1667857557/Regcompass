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
rc_concordance_null_correct <- function(p_rna, p_atac, pool_meta = NULL, stratum_col = NULL) {
  p_rna <- as.matrix(p_rna); p_atac <- as.matrix(p_atac)
  if (!identical(dim(p_rna), dim(p_atac))) stop("`p_rna` and `p_atac` must have identical dimensions.", call. = FALSE)
  concord <- 1 - abs(p_rna - p_atac)
  if (!is.null(stratum_col)) {
    if (is.null(pool_meta) || !"pool_id" %in% colnames(pool_meta) || !stratum_col %in% colnames(pool_meta)) stop("`pool_meta` must contain `pool_id` and `stratum_col`.", call. = FALSE)
    pool_meta <- pool_meta[match(colnames(p_rna), pool_meta$pool_id), , drop = FALSE]
    if (anyNA(pool_meta$pool_id)) stop("`pool_meta` is missing metadata for some matrix columns.", call. = FALSE)
    strata <- as.character(pool_meta[[stratum_col]])
  } else {
    strata <- rep("global", ncol(p_rna))
  }
  out <- concord
  for (st in unique(strata)) {
    cols <- which(strata == st)
    n <- rowSums(is.finite(p_rna[, cols, drop = FALSE]) & is.finite(p_atac[, cols, drop = FALSE]))
    e_null <- 2 / 3 + 1 / (3 * n^2)
    e_null[!is.finite(e_null) | n < 1L] <- NA_real_
    denom <- 1 - e_null
    tmp <- sweep(sweep(concord[, cols, drop = FALSE], 1, e_null, "-"), 1, denom, "/")
    tmp <- pmax(0, pmin(1, tmp))
    zero_power <- !is.finite(denom) | denom <= 0
    if (any(zero_power)) tmp[zero_power, ] <- 0
    tmp[!is.finite(tmp)] <- NA_real_
    out[, cols] <- tmp
  }
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
  out$rel_positive <- pmax(0, out$rho_shrink)
  out$discordance <- abs(pmin(0, out$rho_shrink))
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
  if (is.null(link_conf)) link_conf <- zeros else link_conf <- rc_align_component_matrix(link_conf, base, "link_conf")
  if (is.null(concord_rt_norm)) concord_rt_norm <- zeros else concord_rt_norm <- rc_align_component_matrix(concord_rt_norm, base, "concord_rt_norm")
  if (is.null(rel_rt_pos)) rel_rt_pos <- rep(0, nrow(base))
  if (is.null(gpr_gene_observed)) gpr_gene_observed <- matrix(1, nrow = nrow(base), ncol = ncol(base), dimnames = dimnames(base)) else gpr_gene_observed <- rc_align_component_matrix(gpr_gene_observed, base, "gpr_gene_observed")
  if (is.null(qc)) qc <- rep(1, ncol(base))
  if (!identical(dim(base), dim(link_conf)) || !identical(dim(base), dim(concord_rt_norm)) || !identical(dim(base), dim(as.matrix(gpr_gene_observed)))) {
    stop("All matrix confidence components must have identical dimensions.", call. = FALSE)
  }
  rel_ra_pos <- rc_align_reliability_vector(rel_ra_pos, rownames(base), "rel_ra_pos")
  rel_rt_pos <- rc_align_reliability_vector(rel_rt_pos, rownames(base), "rel_rt_pos")
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

rc_align_reliability_vector <- function(x, genes, name) {
  nms <- names(x)
  x <- as.numeric(x)
  if (!is.null(nms)) {
    names(x) <- nms
    out <- x[match(genes, names(x))]
    out[is.na(out)] <- 0
    return(pmax(0, out))
  }
  if (length(x) != length(genes)) stop("`", name, "` must be named by gene or have one value per gene/row.", call. = FALSE)
  pmax(0, x)
}

rc_align_component_matrix <- function(x, template, name) {
  x <- as.matrix(x)
  if (!is.null(rownames(x)) && !is.null(colnames(x))) {
    if (!all(rownames(template) %in% rownames(x)) || !all(colnames(template) %in% colnames(x))) stop("`", name, "` is missing required genes or pools.", call. = FALSE)
    return(x[rownames(template), colnames(template), drop = FALSE])
  }
  if (!identical(dim(x), dim(template))) stop("`", name, "` must have dimensions matching confidence matrices.", call. = FALSE)
  dimnames(x) <- dimnames(template)
  x
}
