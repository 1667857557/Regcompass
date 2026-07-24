#' Detect a conservative RegCompass worker count
#'
#' Worker discovery honors explicit RegCompass settings before scheduler- or
#' cgroup-aware sources.
#'
#' @param default Fallback worker count.
#' @return A positive integer worker count.
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
  if (!is.finite(cores) || cores < 1L) {
    max(1L, as.integer(default[[1L]]))
  } else {
    max(1L, cores - 1L)
  }
}

.rc_resolve_parallel_backend <- function(
    backend = c("auto", "serial", "snow", "multicore"),
    os_type = .Platform$OS.type) {
  backend <- match.arg(backend)
  os_type <- match.arg(as.character(os_type[[1L]]), c("unix", "windows"))
  if (identical(backend, "auto")) {
    return(if (identical(os_type, "windows")) "snow" else "multicore")
  }
  if (identical(backend, "multicore") && identical(os_type, "windows")) {
    stop(
      "`multicore` is not supported on Windows; use `auto` or `snow`.",
      call. = FALSE
    )
  }
  backend
}

#' Resolve the platform-aware parallel configuration
#'
#' `backend = "auto"` selects a SOCK cluster on Windows and forked multicore
#' workers on Linux/macOS. Sequential execution is used when one worker is
#' requested or BiocParallel is unavailable.
#'
#' @param workers Optional worker count.
#' @param backend Requested backend.
#' @return A list describing requested and resolved execution settings.
#' @export
rc_parallel_config <- function(
    workers = NULL,
    backend = c("auto", "serial", "snow", "multicore")) {
  backend <- match.arg(backend)
  requested_workers <- workers
  if (is.null(workers)) workers <- rc_available_workers(default = 1L)
  workers <- suppressWarnings(as.integer(workers[[1L]]))
  if (!is.finite(workers) || workers < 1L) workers <- 1L
  resolved <- .rc_resolve_parallel_backend(backend)
  available <- requireNamespace("BiocParallel", quietly = TRUE)
  actual <- if (workers < 2L || identical(resolved, "serial") || !available) {
    "serial"
  } else {
    resolved
  }
  list(
    os_type = .Platform$OS.type,
    requested_backend = backend,
    resolved_backend = resolved,
    actual_backend = actual,
    requested_workers = requested_workers,
    workers = if (identical(actual, "serial")) 1L else workers,
    biocparallel_available = available
  )
}

#' Build the default RegCompass parallel backend
#'
#' The backend's task-level progress bar follows
#' `options(RegCompassR.progress = TRUE/FALSE)`.
#'
#' @param workers Optional worker count.
#' @param backend Requested backend.
#' @return A `BiocParallelParam` object or `NULL` for sequential execution.
rc_default_bpparam <- function(
    workers = NULL,
    backend = c("auto", "serial", "snow", "multicore")) {
  config <- rc_parallel_config(workers = workers, backend = backend)
  if (identical(config$actual_backend, "serial")) return(NULL)
  show_progress <- .rc_progress_enabled(
    getOption("RegCompassR.progress", TRUE)
  )

  param <- if (identical(config$actual_backend, "snow")) {
    BiocParallel::SnowParam(
      workers = config$workers,
      type = "SOCK",
      progressbar = show_progress
    )
  } else {
    BiocParallel::MulticoreParam(
      workers = config$workers,
      progressbar = show_progress
    )
  }
  attr(param, "regcompass_parallel_config") <- config
  param
}

#' Apply a function with optional BiocParallel control
#'
#' @param X A vector or list.
#' @param FUN Function applied to each element.
#' @param BPPARAM `NULL`, `FALSE`, or a `BiocParallelParam`.
#' @param ... Additional arguments.
#' @return A list.
rc_parallel_lapply <- function(X, FUN, BPPARAM = NULL, ...) {
  if (identical(BPPARAM, FALSE)) return(lapply(X, FUN, ...))
  if (!is.null(BPPARAM)) {
    if (is.logical(BPPARAM)) {
      stop(
        "`BPPARAM` must be NULL, FALSE, or a BiocParallelParam object; logical TRUE is not valid.",
        call. = FALSE
      )
    }
    if (!requireNamespace("BiocParallel", quietly = TRUE)) {
      stop("BiocParallel must be installed when `BPPARAM` is provided.",
           call. = FALSE)
    }
    if (!methods::is(BPPARAM, "BiocParallelParam")) {
      stop(
        "`BPPARAM` must be NULL, FALSE, or a BiocParallelParam object.",
        call. = FALSE
      )
    }
  }
  if (length(X) <= 1L) return(lapply(X, FUN, ...))
  if (is.null(BPPARAM)) BPPARAM <- rc_default_bpparam()
  if (is.null(BPPARAM)) return(lapply(X, FUN, ...))
  was_started <- isTRUE(BiocParallel::bpisup(BPPARAM))
  if (!was_started) {
    BiocParallel::bpstart(BPPARAM)
    on.exit(BiocParallel::bpstop(BPPARAM), add = TRUE)
  }
  BiocParallel::bplapply(X, FUN, ..., BPPARAM = BPPARAM)
}
