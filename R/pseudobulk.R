#' Sum raw counts by pool
#'
#' This is the only supported expression input for the main Layer 1 workflow:
#' raw counts are summed within each micropool before logCPM normalization.
#' Averaging cell-level residuals, TF-IDF, or imputed expression changes the
#' biological scale and should not be used for reaction capacity potential.
#' @export
rc_pseudobulk_counts <- function(counts, pool_map, fun = "sum", BPPARAM = NULL) {
  if (!identical(fun, "sum")) stop("Main RegCompassR Layer 1 requires `fun = 'sum'`.", call. = FALSE)
  rc_validate_pool_matrix_inputs(counts, pool_map)
  if ("skipped" %in% colnames(pool_map)) pool_map <- pool_map[!pool_map$skipped, , drop = FALSE]
  pool_map <- pool_map[!is.na(pool_map$pool_id), , drop = FALSE]
  pool_ids <- unique(pool_map$pool_id)

  summarize_pool <- function(pid) Matrix::rowSums(counts[, pool_map$cell_id[pool_map$pool_id == pid], drop = FALSE])
  res <- rc_pool_lapply(pool_ids, summarize_pool, BPPARAM = BPPARAM)
  out <- do.call(cbind, res)
  colnames(out) <- pool_ids
  rownames(out) <- rownames(counts)
  out
}

#' Remove zero-library pools before normalization
#' @export
rc_filter_empty_pools <- function(pb_counts, pool_meta) {
  lib <- Matrix::colSums(pb_counts)
  keep <- lib > 0
  if (any(!keep)) warning(sum(!keep), " empty pools removed before normalization", call. = FALSE)
  list(counts = pb_counts[, keep, drop = FALSE],
       pool_meta = pool_meta[match(colnames(pb_counts)[keep], pool_meta$pool_id), , drop = FALSE])
}

#' Pool-level log2(CPM + 1) normalization
#' @export
rc_logcpm <- function(pb_counts, scale_factor = 1e6) {
  lib <- Matrix::colSums(pb_counts)
  if (any(lib <= 0)) stop("Empty pools detected. Run rc_filter_empty_pools() first.", call. = FALSE)
  norm <- t(t(pb_counts) / lib) * scale_factor
  log1p(norm) / log(2)
}

#' Build one-row-per-pool metadata
#' @export
rc_build_pool_metadata <- function(pool_map, meta = NULL) {
  keep <- !is.na(pool_map$pool_id)
  if ("skipped" %in% colnames(pool_map)) keep <- keep & !pool_map$skipped
  x <- pool_map[keep, , drop = FALSE]
  cols <- setdiff(colnames(x), "cell_id")
  out <- x[!duplicated(x$pool_id), cols, drop = FALSE]
  out$n_cells <- as.integer(tabulate(match(x$pool_id, out$pool_id), nbins = nrow(out)))
  rownames(out) <- NULL
  out
}

#' Pool-level RNA detection rates for confidence only
#' @export
rc_pool_detection <- function(counts, pool_map, BPPARAM = NULL) {
  rc_validate_pool_matrix_inputs(counts, pool_map)
  if ("skipped" %in% colnames(pool_map)) pool_map <- pool_map[!pool_map$skipped, , drop = FALSE]
  pool_map <- pool_map[!is.na(pool_map$pool_id), , drop = FALSE]
  pool_ids <- unique(pool_map$pool_id)
  summarize_pool <- function(pid) Matrix::rowMeans(counts[, pool_map$cell_id[pool_map$pool_id == pid], drop = FALSE] > 0)
  res <- rc_pool_lapply(pool_ids, summarize_pool, BPPARAM = BPPARAM)
  out <- do.call(cbind, res)
  colnames(out) <- pool_ids
  rownames(out) <- rownames(counts)
  out
}

rc_validate_pool_matrix_inputs <- function(mat, pool_map) {
  if (is.null(dim(mat)) || length(dim(mat)) != 2L) stop("`counts` must be a two-dimensional feature-by-cell matrix.", call. = FALSE)
  if (is.null(colnames(mat)) || anyNA(colnames(mat)) || any(!nzchar(colnames(mat)))) stop("Input matrix must have non-empty cell IDs in colnames().", call. = FALSE)
  if (!is.data.frame(pool_map)) stop("`pool_map` must be a data.frame.", call. = FALSE)
  missing_pool_cols <- setdiff(c("pool_id", "cell_id"), colnames(pool_map))
  if (length(missing_pool_cols) > 0L) stop("`pool_map` is missing required columns: ", paste(missing_pool_cols, collapse = ", "), call. = FALSE)
  active <- pool_map
  if ("skipped" %in% colnames(active)) active <- active[!active$skipped, , drop = FALSE]
  active <- active[!is.na(active$pool_id), , drop = FALSE]
  if (anyNA(active$cell_id)) stop("`pool_map$cell_id` must not contain NA values.", call. = FALSE)
  missing_cells <- setdiff(active$cell_id, colnames(mat))
  if (length(missing_cells) > 0L) stop("Some pool_map cell IDs are absent from matrix columns: ", paste(utils::head(missing_cells), collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

rc_pool_lapply <- function(X, FUN, BPPARAM = NULL) {
  if (!is.null(BPPARAM)) {
    if (!requireNamespace("BiocParallel", quietly = TRUE)) stop("BiocParallel must be installed when `BPPARAM` is provided.", call. = FALSE)
    return(BiocParallel::bplapply(X, FUN, BPPARAM = BPPARAM))
  }
  lapply(X, FUN)
}
