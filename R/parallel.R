#' Detect a conservative RegCompass worker count
#'
#' Worker discovery honors explicit RegCompass settings before scheduler- or
#' cgroup-aware sources. This avoids over-subscribing containers where
#' `parallel::detectCores()` reports host CPUs instead of the current process
#' allocation.
#'
#' @param default Fallback worker count when no source can be detected.
#' @return A positive integer worker count.
#' @export
rc_available_workers <- function(default = 1L) {
  vals <- c(
    getOption("RegCompassR.workers", NA),
    Sys.getenv("REGCOMPASS_WORKERS", unset = NA_character_),
    Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA_character_),
    Sys.getenv("NSLOTS", unset = NA_character_)
  )
  vals <- suppressWarnings(as.integer(vals))
  vals <- vals[is.finite(vals) & vals >= 1L]
  if (length(vals)) return(max(1L, vals[[1L]]))

  if (requireNamespace("future", quietly = TRUE)) {
    fc <- tryCatch(future::availableCores(), error = function(e) NA_integer_)
    fc <- suppressWarnings(as.integer(fc[[1L]]))
    if (is.finite(fc) && fc >= 1L) return(fc)
  }

  cores <- parallel::detectCores(logical = TRUE)
  cores <- suppressWarnings(as.integer(cores[[1L]]))
  if (!is.finite(cores) || cores < 1L) max(1L, as.integer(default[[1L]])) else max(1L, cores - 1L)
}

#' Build the default RegCompass parallel backend
#'
#' Expensive pool-, reaction-, and bootstrap-level loops use this helper when
#' callers do not provide an explicit `BiocParallelParam`. By default it uses a
#' conservative worker count from `rc_available_workers()`. Set
#' `options(RegCompassR.workers = 1)` or `REGCOMPASS_WORKERS=1` to force
#' sequential execution.
#'
#' @param workers Optional worker count. Defaults to option/env/autodetection.
#' @return A `BiocParallelParam` object when BiocParallel is installed and more
#' than one worker is requested; otherwise `NULL` for sequential execution.
#' @export
rc_default_bpparam <- function(workers = NULL) {
  if (is.null(workers)) workers <- rc_available_workers(default = 1L)
  workers <- suppressWarnings(as.integer(workers[[1L]]))
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
