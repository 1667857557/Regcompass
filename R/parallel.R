#' Build the default RegCompass parallel backend
#'
#' Expensive pool-, reaction-, and bootstrap-level loops use this helper when
#' callers do not provide an explicit `BiocParallelParam`. By default it uses all
#' available cores except one, capped by `options(RegCompassR.workers = ...)` or
#' the `REGCOMPASS_WORKERS` environment variable. Set either value to `1` to force
#' sequential execution.
#'
#' @param workers Optional worker count. Defaults to option/env autodetection.
#' @return A `BiocParallelParam` object when BiocParallel is installed and more
#' than one worker is requested; otherwise `NULL` for sequential execution.
#' @export
rc_default_bpparam <- function(workers = NULL) {
  if (is.null(workers)) {
    opt_workers <- getOption("RegCompassR.workers", NULL)
    env_workers <- Sys.getenv("REGCOMPASS_WORKERS", unset = NA_character_)
    if (!is.null(opt_workers)) {
      workers <- opt_workers
    } else if (!is.na(env_workers) && nzchar(env_workers)) {
      workers <- env_workers
    } else {
      cores <- parallel::detectCores(logical = TRUE)
      workers <- if (is.na(cores)) 1L else max(1L, cores - 1L)
    }
  }
  workers <- suppressWarnings(as.integer(workers[[1]]))
  if (is.na(workers) || workers < 2L) return(NULL)
  if (!requireNamespace("BiocParallel", quietly = TRUE)) return(NULL)

  if (.Platform$OS.type == "windows") {
    BiocParallel::SnowParam(workers = workers, type = "SOCK")
  } else {
    BiocParallel::MulticoreParam(workers = workers)
  }
}

#' Apply a function with optional BiocParallel control
#'
#' @param X A vector or list to iterate over.
#' @param FUN Function applied to each element of `X`.
#' @param BPPARAM Optional `BiocParallelParam`. If `NULL`, RegCompass attempts to
#' use `rc_default_bpparam()` for multi-core work and falls back to base
#' `lapply()` when BiocParallel is unavailable or only one worker is requested.
#' Pass `BPPARAM = FALSE` to force sequential execution.
#' @param ... Additional arguments passed to `FUN`.
#'
#' @return A list with one element per `X`.
#' @export
rc_parallel_lapply <- function(X, FUN, BPPARAM = NULL, ...) {
  if (identical(BPPARAM, FALSE) || length(X) <= 1L) {
    return(lapply(X, FUN, ...))
  }
  if (is.null(BPPARAM)) BPPARAM <- rc_default_bpparam()
  if (is.null(BPPARAM)) return(lapply(X, FUN, ...))
  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    stop("BiocParallel must be installed when `BPPARAM` is provided.", call. = FALSE)
  }
  BiocParallel::bplapply(X, FUN, ..., BPPARAM = BPPARAM)
}
