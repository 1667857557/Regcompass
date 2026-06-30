#' Compute pool-level sparse means
#'
#' `rc_pool_mean()` summarizes a gene-by-cell matrix into a gene-by-pool matrix
#' using sparse-aware row means. Use normalized data or residual matrices for
#' expression scores; raw counts should be reserved for detection-rate summaries.
#'
#' @param mat A feature-by-cell matrix, typically sparse, with cell IDs in
#' `colnames(mat)`.
#' @param pool_map A pool assignment data.frame containing `pool_id` and
#' `cell_id`, usually returned by [rc_make_pools()].
#' @param BPPARAM Optional `BiocParallelParam` object. When provided and
#' BiocParallel is installed, pools are summarized with
#' `BiocParallel::bplapply()`; otherwise base `lapply()` is used.
#'
#' @return A numeric feature-by-pool matrix with pool IDs as column names.
#' @export
rc_pool_mean <- function(mat, pool_map, BPPARAM = NULL) {
  rc_validate_pool_matrix_inputs(mat, pool_map)
  pool_ids <- unique(pool_map$pool_id)

  summarize_pool <- function(pid) {
    cells <- pool_map$cell_id[pool_map$pool_id == pid]
    Matrix::rowMeans(mat[, cells, drop = FALSE])
  }

  res <- rc_pool_lapply(pool_ids, summarize_pool, BPPARAM = BPPARAM)
  out <- do.call(cbind, res)
  colnames(out) <- pool_ids
  rownames(out) <- rownames(mat)
  out
}

#' Compute pool-level detection rates
#'
#' `rc_pool_detection()` summarizes raw counts into the fraction of cells in each
#' pool where a feature is detected (`counts > 0`). This detection rate is the
#' downstream dropout-aware quantity `q_{g,p}` in the development specification.
#'
#' @param counts A raw feature-by-cell count matrix, typically sparse, with cell
#' IDs in `colnames(counts)`.
#' @inheritParams rc_pool_mean
#'
#' @return A numeric feature-by-pool matrix of detection fractions in `[0, 1]`.
#' @export
rc_pool_detection <- function(counts, pool_map, BPPARAM = NULL) {
  rc_validate_pool_matrix_inputs(counts, pool_map)
  pool_ids <- unique(pool_map$pool_id)

  summarize_pool <- function(pid) {
    cells <- pool_map$cell_id[pool_map$pool_id == pid]
    Matrix::rowMeans(counts[, cells, drop = FALSE] > 0)
  }

  res <- rc_pool_lapply(pool_ids, summarize_pool, BPPARAM = BPPARAM)
  out <- do.call(cbind, res)
  colnames(out) <- pool_ids
  rownames(out) <- rownames(counts)
  out
}

rc_validate_pool_matrix_inputs <- function(mat, pool_map) {
  if (is.null(dim(mat)) || length(dim(mat)) != 2L) {
    stop("`mat`/`counts` must be a two-dimensional feature-by-cell matrix.", call. = FALSE)
  }
  if (is.null(colnames(mat)) || anyNA(colnames(mat)) || any(!nzchar(colnames(mat)))) {
    stop("Input matrix must have non-empty cell IDs in colnames().", call. = FALSE)
  }
  if (!is.data.frame(pool_map)) {
    stop("`pool_map` must be a data.frame.", call. = FALSE)
  }
  missing_pool_cols <- setdiff(c("pool_id", "cell_id"), colnames(pool_map))
  if (length(missing_pool_cols) > 0) {
    stop("`pool_map` is missing required columns: ", paste(missing_pool_cols, collapse = ", "), call. = FALSE)
  }
  if (anyNA(pool_map$pool_id) || anyNA(pool_map$cell_id)) {
    stop("`pool_map$pool_id` and `pool_map$cell_id` must not contain NA values.", call. = FALSE)
  }
  missing_cells <- setdiff(pool_map$cell_id, colnames(mat))
  if (length(missing_cells) > 0) {
    stop("Some pool_map cell IDs are absent from matrix columns: ", paste(utils::head(missing_cells), collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

rc_pool_lapply <- function(X, FUN, BPPARAM = NULL) {
  if (!is.null(BPPARAM)) {
    if (!requireNamespace("BiocParallel", quietly = TRUE)) {
      stop("BiocParallel must be installed when `BPPARAM` is provided.", call. = FALSE)
    }
    return(BiocParallel::bplapply(X, FUN, BPPARAM = BPPARAM))
  }
  lapply(X, FUN)
}
