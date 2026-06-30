#' Apply a function with optional BiocParallel control
#'
#' @param X A vector or list to iterate over.
#' @param FUN Function applied to each element of `X`.
#' @param BPPARAM Optional `BiocParallelParam`. If `NULL`, base `lapply()` is
#' used; otherwise `BiocParallel::bplapply()` is used.
#' @param ... Additional arguments passed to `FUN`.
#'
#' @return A list with one element per `X`.
#' @export
rc_parallel_lapply <- function(X, FUN, BPPARAM = NULL, ...) {
  if (is.null(BPPARAM)) {
    return(lapply(X, FUN, ...))
  }
  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    stop("BiocParallel must be installed when `BPPARAM` is provided.", call. = FALSE)
  }
  BiocParallel::bplapply(X, FUN, ..., BPPARAM = BPPARAM)
}
