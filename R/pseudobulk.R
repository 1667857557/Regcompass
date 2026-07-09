#' Pseudobulk raw counts by pool
rc_unit_bulk_counts <- function(counts, unit_map, fun = c("sum", "mean"), BPPARAM = NULL) {
  fun <- match.arg(fun)
  rc_validate_unit_matrix_inputs(counts, unit_map)
  unit_map <- rc_filter_active_unit_map(unit_map)
  pool_ids <- unique(unit_map$pool_id)
  summarize_pool <- function(pid) {
    cells <- unit_map$cell_id[unit_map$pool_id == pid]
    x <- counts[, cells, drop = FALSE]
    if (fun == "sum") Matrix::rowSums(x) else Matrix::rowMeans(x)
  }
  res <- rc_internal_lapply(pool_ids, summarize_pool, BPPARAM = BPPARAM)
  out <- do.call(cbind, res)
  colnames(out) <- pool_ids; rownames(out) <- rownames(counts)
  out
}

#' Filter empty pseudobulk pools before normalization
rc_filter_empty_units <- function(pb_counts, unit_meta) {
  lib <- Matrix::colSums(pb_counts)
  keep <- lib > 0
  if (any(!keep)) warning(sum(!keep), " empty units removed before normalization", call. = FALSE)
  list(counts = pb_counts[, keep, drop = FALSE],
       unit_meta = unit_meta[match(colnames(pb_counts)[keep], unit_meta$pool_id), , drop = FALSE])
}

#' Meta-cell-level log2(CPM + 1) normalization
#' @export
rc_logcpm <- function(pb_counts, scale_factor = 1e6) {
  pb_counts <- .rc_as_sparse(pb_counts)
  lib <- Matrix::colSums(pb_counts)
  if (any(lib <= 0)) stop("Empty units detected. Run rc_filter_empty_units() first.", call. = FALSE)
  norm <- Matrix::t(Matrix::t(pb_counts) / lib) * scale_factor
  log1p(norm) / log(2)
}

#' Build one-row-per-pool metadata
rc_build_unit_metadata <- function(unit_map, meta = NULL) {
  x <- rc_filter_active_unit_map(unit_map)
  cols <- setdiff(colnames(x), "cell_id")
  out <- x[!duplicated(x$pool_id), cols, drop = FALSE]
  out$n_cells <- as.integer(tabulate(match(x$pool_id, out$pool_id), nbins = nrow(out)))
  rownames(out) <- NULL
  out
}

#' Compute meta-cell-level sparse means
rc_unit_mean <- function(mat, unit_map, BPPARAM = NULL) {
  rc_validate_unit_matrix_inputs(mat, unit_map)
  unit_map <- rc_filter_active_unit_map(unit_map)
  pool_ids <- unique(unit_map$pool_id)
  summarize_pool <- function(pid) Matrix::rowMeans(mat[, unit_map$cell_id[unit_map$pool_id == pid], drop = FALSE])
  res <- rc_internal_lapply(pool_ids, summarize_pool, BPPARAM = BPPARAM)
  out <- do.call(cbind, res); colnames(out) <- pool_ids; rownames(out) <- rownames(mat); out
}

#' Compute meta-cell-level detection rates
rc_unit_detection <- function(counts, unit_map, BPPARAM = NULL) {
  rc_validate_unit_matrix_inputs(counts, unit_map)
  unit_map <- rc_filter_active_unit_map(unit_map)
  pool_ids <- unique(unit_map$pool_id)
  summarize_pool <- function(pid) Matrix::rowMeans(counts[, unit_map$cell_id[unit_map$pool_id == pid], drop = FALSE] > 0)
  res <- rc_internal_lapply(pool_ids, summarize_pool, BPPARAM = BPPARAM)
  out <- do.call(cbind, res); colnames(out) <- pool_ids; rownames(out) <- rownames(counts); out
}

rc_validate_unit_matrix_inputs <- function(mat, unit_map) {
  if (is.null(dim(mat)) || length(dim(mat)) != 2L) stop("`mat`/`counts` must be a two-dimensional feature-by-cell matrix.", call. = FALSE)
  if (is.null(colnames(mat)) || anyNA(colnames(mat)) || any(!nzchar(colnames(mat)))) stop("Input matrix must have non-empty cell IDs in colnames().", call. = FALSE)
  if (!is.data.frame(unit_map)) stop("`unit_map` must be a data.frame.", call. = FALSE)
  missing_pool_cols <- setdiff(c("pool_id", "cell_id"), colnames(unit_map))
  if (length(missing_pool_cols) > 0) stop("`unit_map` is missing required columns: ", paste(missing_pool_cols, collapse = ", "), call. = FALSE)
  active <- rc_filter_active_unit_map(unit_map)
  if (anyNA(active$cell_id)) stop("`unit_map$cell_id` must not contain NA values.", call. = FALSE)
  missing_cells <- setdiff(active$cell_id, colnames(mat))
  if (length(missing_cells) > 0L) stop("Some unit_map cell IDs are absent from matrix columns: ", paste(utils::head(missing_cells), collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

rc_filter_active_unit_map <- function(unit_map) {
  keep <- !is.na(unit_map$pool_id)
  if ("skipped" %in% colnames(unit_map)) keep <- keep & !(unit_map$skipped %in% TRUE)
  unit_map[keep, , drop = FALSE]
}

rc_internal_lapply <- function(X, FUN, BPPARAM = NULL) {
  rc_parallel_lapply(X, FUN, BPPARAM = BPPARAM)
}

#' Filter ATAC peaks detected in at least min_pools and compute pool logCPM
rc_atac_unit_logcpm <- function(atac_counts, unit_map, min_pools = 3, BPPARAM = NULL) {
  pb <- rc_unit_bulk_counts(atac_counts, unit_map, fun = "sum", BPPARAM = BPPARAM)
  detected <- Matrix::rowSums(pb > 0) >= min_pools
  pb <- pb[detected, , drop = FALSE]
  filtered <- rc_filter_empty_units(pb, rc_build_unit_metadata(unit_map))
  rc_logcpm(filtered$counts)
}

#' Check pseudobulk columns against manual pool sums
#'
#' This lightweight runtime sanity check verifies that selected pseudobulk
#' columns equal manual row sums over the active pool-map cells.
rc_check_unit_mapping <- function(counts, unit_map, pb, n_check = 5) {
  rc_validate_unit_matrix_inputs(counts, unit_map)
  active <- rc_filter_active_unit_map(unit_map)
  pids <- unique(active$pool_id)
  pids <- utils::head(pids, n_check)
  missing_pools <- setdiff(pids, colnames(pb))
  if (length(missing_pools) > 0L) stop("`pb` is missing unit columns: ", paste(missing_pools, collapse = ", "), call. = FALSE)
  for (pid in pids) {
    cells <- active$cell_id[active$pool_id == pid]
    manual <- Matrix::rowSums(counts[, cells, drop = FALSE])
    if (!isTRUE(all.equal(as.numeric(pb[, pid]), as.numeric(manual)))) {
      stop("Pseudobulk mapping check failed for pool: ", pid, call. = FALSE)
    }
  }
  TRUE
}
