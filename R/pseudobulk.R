#' Pseudobulk raw counts by pool
#' @export
rc_pseudobulk_counts <- function(counts, pool_map, fun = c("sum", "mean"), BPPARAM = NULL) {
  fun <- match.arg(fun)
  rc_validate_pool_matrix_inputs(counts, pool_map)
  pool_map <- rc_filter_active_pool_map(pool_map)
  pool_ids <- unique(pool_map$pool_id)
  summarize_pool <- function(pid) {
    cells <- pool_map$cell_id[pool_map$pool_id == pid]
    x <- counts[, cells, drop = FALSE]
    if (fun == "sum") Matrix::rowSums(x) else Matrix::rowMeans(x)
  }
  res <- rc_pool_lapply(pool_ids, summarize_pool, BPPARAM = BPPARAM)
  out <- do.call(cbind, res)
  colnames(out) <- pool_ids; rownames(out) <- rownames(counts)
  out
}

#' Filter empty pseudobulk pools before normalization
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
  x <- rc_filter_active_pool_map(pool_map)
  cols <- setdiff(colnames(x), "cell_id")
  out <- x[!duplicated(x$pool_id), cols, drop = FALSE]
  out$n_cells <- as.integer(tabulate(match(x$pool_id, out$pool_id), nbins = nrow(out)))
  rownames(out) <- NULL
  out
}

#' Compute pool-level sparse means
#' @export
rc_pool_mean <- function(mat, pool_map, BPPARAM = NULL) {
  rc_validate_pool_matrix_inputs(mat, pool_map)
  pool_map <- rc_filter_active_pool_map(pool_map)
  pool_ids <- unique(pool_map$pool_id)
  summarize_pool <- function(pid) Matrix::rowMeans(mat[, pool_map$cell_id[pool_map$pool_id == pid], drop = FALSE])
  res <- rc_pool_lapply(pool_ids, summarize_pool, BPPARAM = BPPARAM)
  out <- do.call(cbind, res); colnames(out) <- pool_ids; rownames(out) <- rownames(mat); out
}

#' Compute pool-level detection rates
#' @export
rc_pool_detection <- function(counts, pool_map, BPPARAM = NULL) {
  rc_validate_pool_matrix_inputs(counts, pool_map)
  pool_map <- rc_filter_active_pool_map(pool_map)
  pool_ids <- unique(pool_map$pool_id)
  summarize_pool <- function(pid) Matrix::rowMeans(counts[, pool_map$cell_id[pool_map$pool_id == pid], drop = FALSE] > 0)
  res <- rc_pool_lapply(pool_ids, summarize_pool, BPPARAM = BPPARAM)
  out <- do.call(cbind, res); colnames(out) <- pool_ids; rownames(out) <- rownames(counts); out
}

rc_validate_pool_matrix_inputs <- function(mat, pool_map) {
  if (is.null(dim(mat)) || length(dim(mat)) != 2L) stop("`mat`/`counts` must be a two-dimensional feature-by-cell matrix.", call. = FALSE)
  if (is.null(colnames(mat)) || anyNA(colnames(mat)) || any(!nzchar(colnames(mat)))) stop("Input matrix must have non-empty cell IDs in colnames().", call. = FALSE)
  if (!is.data.frame(pool_map)) stop("`pool_map` must be a data.frame.", call. = FALSE)
  missing_pool_cols <- setdiff(c("pool_id", "cell_id"), colnames(pool_map))
  if (length(missing_pool_cols) > 0) stop("`pool_map` is missing required columns: ", paste(missing_pool_cols, collapse = ", "), call. = FALSE)
  active <- rc_filter_active_pool_map(pool_map)
  if (anyNA(active$cell_id)) stop("`pool_map$cell_id` must not contain NA values.", call. = FALSE)
  missing_cells <- setdiff(active$cell_id, colnames(mat))
  if (length(missing_cells) > 0L) stop("Some pool_map cell IDs are absent from matrix columns: ", paste(utils::head(missing_cells), collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

rc_filter_active_pool_map <- function(pool_map) {
  keep <- !is.na(pool_map$pool_id)
  if ("skipped" %in% colnames(pool_map)) keep <- keep & !(pool_map$skipped %in% TRUE)
  pool_map[keep, , drop = FALSE]
}

rc_pool_lapply <- function(X, FUN, BPPARAM = NULL) {
  rc_parallel_lapply(X, FUN, BPPARAM = BPPARAM)
}

#' Filter ATAC peaks detected in at least min_pools and compute pool logCPM
#' @export
rc_atac_pool_logcpm <- function(atac_counts, pool_map, min_pools = 3, BPPARAM = NULL) {
  pb <- rc_pseudobulk_counts(atac_counts, pool_map, fun = "sum", BPPARAM = BPPARAM)
  detected <- Matrix::rowSums(pb > 0) >= min_pools
  pb <- pb[detected, , drop = FALSE]
  filtered <- rc_filter_empty_pools(pb, rc_build_pool_metadata(pool_map))
  rc_logcpm(filtered$counts)
}

#' Check pseudobulk columns against manual pool sums
#'
#' This lightweight runtime sanity check verifies that selected pseudobulk
#' columns equal manual row sums over the active pool-map cells.
#' @export
rc_check_pseudobulk_mapping <- function(counts, pool_map, pb, n_check = 5) {
  rc_validate_pool_matrix_inputs(counts, pool_map)
  active <- rc_filter_active_pool_map(pool_map)
  pids <- unique(active$pool_id)
  pids <- utils::head(pids, n_check)
  missing_pools <- setdiff(pids, colnames(pb))
  if (length(missing_pools) > 0L) stop("`pb` is missing pool columns: ", paste(missing_pools, collapse = ", "), call. = FALSE)
  for (pid in pids) {
    cells <- active$cell_id[active$pool_id == pid]
    manual <- Matrix::rowSums(counts[, cells, drop = FALSE])
    if (!isTRUE(all.equal(as.numeric(pb[, pid]), as.numeric(manual)))) {
      stop("Pseudobulk mapping check failed for pool: ", pid, call. = FALSE)
    }
  }
  TRUE
}
