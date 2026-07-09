#' Stratum-wise empirical percentiles for feature-by-pool matrices
#'
#' Percentiles are computed within each requested unit stratum so cell-type or
#' condition composition does not dominate RNA/ATAC confidence components.
#' @export
rc_percentile_by_stratum <- function(x, unit_meta = NULL, stratum_col = NULL) {
  x <- as.matrix(x)
  if (is.null(colnames(x))) stop("`x` must have unit IDs in colnames().", call. = FALSE)
  if (is.null(stratum_col)) {
    strata <- rep("global", ncol(x))
  } else {
    if (is.null(unit_meta) || !"pool_id" %in% colnames(unit_meta) || !stratum_col %in% colnames(unit_meta)) {
      stop("`unit_meta` must contain `pool_id` and `stratum_col`.", call. = FALSE)
    }
    unit_meta <- unit_meta[match(colnames(x), unit_meta$pool_id), , drop = FALSE]
    if (anyNA(unit_meta$pool_id)) stop("`unit_meta` is missing metadata for some matrix columns.", call. = FALSE)
    strata <- as.character(unit_meta[[stratum_col]])
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
    out[ok] <- NA_real_
    return(out)
  }
  out[ok] <- (rank(v[ok], ties.method = "average") - 1) / (n - 1)
  out
}

#' RNA-ATAC concordance with discrete-rank null correction
#' @export
rc_concordance_null_correct <- function(p_rna, p_atac, unit_meta = NULL, stratum_col = NULL) {
  p_rna <- as.matrix(p_rna); p_atac <- as.matrix(p_atac)
  if (!identical(dim(p_rna), dim(p_atac))) stop("`p_rna` and `p_atac` must have identical dimensions.", call. = FALSE)
  concord <- 1 - abs(p_rna - p_atac)
  if (!is.null(stratum_col)) {
    if (is.null(unit_meta) || !"pool_id" %in% colnames(unit_meta) || !stratum_col %in% colnames(unit_meta)) stop("`unit_meta` must contain `pool_id` and `stratum_col`.", call. = FALSE)
    unit_meta <- unit_meta[match(colnames(p_rna), unit_meta$pool_id), , drop = FALSE]
    if (anyNA(unit_meta$pool_id)) stop("`unit_meta` is missing metadata for some matrix columns.", call. = FALSE)
    strata <- as.character(unit_meta[[stratum_col]])
  } else {
    strata <- rep("global", ncol(p_rna))
  }
  out <- concord
  for (st in unique(strata)) {
    cols <- which(strata == st)
    n <- rowSums(is.finite(p_rna[, cols, drop = FALSE]) & is.finite(p_atac[, cols, drop = FALSE]))
    e_null <- 2 / 3 - 1 / (3 * n)
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
                               rel_rt_pos = NULL,
                               return_components = FALSE) {
  base <- as.matrix(concord_ra_norm)
  det_rna <- as.matrix(det_rna)
  if (!identical(dim(base), dim(det_rna))) stop("`concord_ra_norm` and `det_rna` must have identical dimensions.", call. = FALSE)
  missing_components <- character(0)
  if (is.null(link_conf)) {
    missing_components <- c(missing_components, "link_conf")
  } else {
    link_conf <- rc_align_component_matrix(link_conf, base, "link_conf")
  }
  if (is.null(concord_rt_norm)) {
    missing_components <- c(missing_components, "concord_rt_norm")
  } else {
    concord_rt_norm <- rc_align_component_matrix(concord_rt_norm, base, "concord_rt_norm")
  }
  if (is.null(rel_rt_pos)) {
    missing_components <- c(missing_components, "rel_rt_pos")
  }
  if (is.null(gpr_gene_observed)) {
    missing_components <- c(missing_components, "gpr_gene_observed")
  } else {
    gpr_gene_observed <- rc_align_component_matrix(gpr_gene_observed, base, "gpr_gene_observed")
  }
  if (is.null(qc)) {
    missing_components <- c(missing_components, "qc")
  }
  rel_ra_pos <- rc_align_reliability_vector(rel_ra_pos, rownames(base), "rel_ra_pos")
  if (!is.null(rel_rt_pos)) rel_rt_pos <- rc_align_reliability_vector(rel_rt_pos, rownames(base), "rel_rt_pos")

  components <- list(
    ra = list(weight = 0.25, value = sweep(rc_clamp01_matrix(base), 1, rel_ra_pos, "*")),
    det = list(weight = 0.20, value = rc_clamp01_matrix(det_rna))
  )
  if (!is.null(concord_rt_norm) && !is.null(rel_rt_pos)) {
    components$rt <- list(weight = 0.15, value = sweep(rc_clamp01_matrix(concord_rt_norm), 1, rel_rt_pos, "*"))
  }
  if (!is.null(link_conf)) {
    components$link <- list(weight = 0.15, value = rc_clamp01_matrix(link_conf))
  }
  if (!is.null(qc)) {
    if (length(qc) != ncol(base)) stop("`qc` must have one value per pool/column.", call. = FALSE)
    qc_mat <- matrix(pmax(0, pmin(1, as.numeric(qc))), nrow = nrow(base), ncol = ncol(base), byrow = TRUE)
    components$qc <- list(weight = 0.15, value = qc_mat)
  }
  if (!is.null(gpr_gene_observed)) {
    components$gpr_observed <- list(weight = 0.10, value = rc_clamp01_matrix(gpr_gene_observed))
  }

  num <- matrix(0, nrow = nrow(base), ncol = ncol(base), dimnames = dimnames(base))
  den <- matrix(0, nrow = nrow(base), ncol = ncol(base), dimnames = dimnames(base))
  for (component in components) {
    val <- component$value
    if (!identical(dim(base), dim(val))) stop("All matrix confidence components must have identical dimensions.", call. = FALSE)
    ok <- is.finite(val)
    num[ok] <- num[ok] + component$weight * val[ok]
    den[ok] <- den[ok] + component$weight
  }
  conf <- num / den
  conf[den == 0] <- NA_real_
  conf <- rc_clamp01_matrix(conf)
  attr(conf, "confidence_component_missing_flag") <- length(missing_components) > 0L
  attr(conf, "missing_components") <- missing_components
  if (isTRUE(return_components)) {
    component_values <- lapply(components, `[[`, "value")
    return(list(
      gene_confidence = conf,
      confidence = conf,
      ra_component = component_values$ra,
      det_component = component_values$det,
      link_component = if (!is.null(component_values$link)) component_values$link else NULL,
      qc_component = if (!is.null(component_values$qc)) component_values$qc else NULL,
      gpr_observed_component = if (!is.null(component_values$gpr_observed)) component_values$gpr_observed else NULL,
      rel_ra_pos = rel_ra_pos,
      concord_ra_norm = rc_clamp01_matrix(base),
      components = component_values,
      component_weights = vapply(components, `[[`, numeric(1), "weight"),
      missing_components = missing_components
    ))
  }
  conf
}

rc_clamp01_matrix <- function(x) {
  x <- as.matrix(x)
  x[] <- pmax(0, pmin(1, x))
  x
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

#' Stratum-aware ATAC link confidence from normalized peak accessibility percentiles
#' @export
rc_link_confidence_by_stratum <- function(p_atac_peak,
                                          peak_gene_links,
                                          unit_meta,
                                          link_stratum_cols = "cell_type") {
  if (!"link_stratum" %in% colnames(peak_gene_links)) stop("`peak_gene_links` must contain `link_stratum`.", call. = FALSE)
  missing_cols <- setdiff(link_stratum_cols, colnames(unit_meta))
  if (length(missing_cols) > 0L) stop("`unit_meta` missing link stratum columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  unit_stratum <- interaction(unit_meta[, link_stratum_cols, drop = FALSE], sep = "|", drop = TRUE)
  names(unit_stratum) <- as.character(unit_meta$pool_id)
  out <- matrix(NA_real_, nrow = length(unique(peak_gene_links$gene)), ncol = ncol(p_atac_peak), dimnames = list(sort(unique(peak_gene_links$gene)), colnames(p_atac_peak)))
  for (st in unique(as.character(peak_gene_links$link_stratum))) {
    cols <- names(unit_stratum)[as.character(unit_stratum) == st]
    cols <- intersect(cols, colnames(p_atac_peak))
    if (length(cols) == 0L) next
    links_st <- peak_gene_links[peak_gene_links$link_stratum == st, , drop = FALSE]
    conf_st <- rc_link_confidence(p_atac_peak[, cols, drop = FALSE], links_st)
    common_genes <- intersect(rownames(out), rownames(conf_st))
    out[common_genes, cols] <- conf_st[common_genes, cols, drop = FALSE]
  }
  out
}
