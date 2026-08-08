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

.rc_validate_worker_limit <- function(workers = NULL, argument = "workers") {
  if (is.null(workers)) return(rc_available_workers(default = 1L))
  if (length(workers) != 1L || is.na(workers) || !is.finite(workers)) {
    stop("`", argument, "` must be one positive integer.", call. = FALSE)
  }
  workers <- suppressWarnings(as.integer(workers))
  if (workers < 1L) {
    stop("`", argument, "` must be at least 1.", call. = FALSE)
  }
  workers
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

.rc_bpparam_backend <- function(BPPARAM) {
  if (identical(BPPARAM, FALSE) || is.null(BPPARAM)) return("serial")
  if (!requireNamespace("BiocParallel", quietly = TRUE) ||
      !methods::is(BPPARAM, "BiocParallelParam")) return("auto")
  if (methods::is(BPPARAM, "SnowParam")) return("snow")
  if (methods::is(BPPARAM, "MulticoreParam")) return("multicore")
  "auto"
}

.rc_bpparam_worker_limit <- function(BPPARAM, default = 1L) {
  if (identical(BPPARAM, FALSE) || is.null(BPPARAM)) {
    return(max(1L, as.integer(default)))
  }
  recorded <- suppressWarnings(as.integer(
    attr(BPPARAM, "regcompass_worker_limit") %||% NA_integer_
  ))
  if (is.finite(recorded) && recorded >= 1L) return(recorded)
  if (requireNamespace("BiocParallel", quietly = TRUE) &&
      methods::is(BPPARAM, "BiocParallelParam")) {
    return(max(1L, as.integer(BiocParallel::bpnworkers(BPPARAM))))
  }
  max(1L, as.integer(default))
}

#' Resolve the platform-aware parallel configuration
#'
#' `backend = "auto"` selects a SOCK cluster on Windows and forked multicore
#' workers on Linux/macOS. `workers` is the single RegCompass-wide worker cap:
#' individual dispatches may use fewer workers when fewer independent tasks are
#' available, but no package-managed dispatch may exceed this value.
#'
#' @param workers Optional total worker cap. When `NULL`, RegCompass detects the
#'   available worker count once from scheduler/cgroup/local resources.
#' @param backend Requested backend.
#' @return A list describing requested and resolved execution settings.
#' @export
rc_parallel_config <- function(
    workers = NULL,
    backend = c("auto", "serial", "snow", "multicore")) {
  backend <- match.arg(backend)
  requested_workers <- workers
  workers <- .rc_validate_worker_limit(workers, argument = "workers")
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
    worker_limit = workers,
    workers = if (identical(actual, "serial")) 1L else workers,
    biocparallel_available = available
  )
}

#' Build the default RegCompass parallel backend
#'
#' The backend's task-level progress bar follows
#' `options(RegCompassR.progress = TRUE/FALSE)`. The returned parameter records
#' the total RegCompass worker cap. Dispatchers create smaller short-lived pools
#' when the current task count is below that cap.
#'
#' @param workers Optional total worker cap.
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
  attr(param, "regcompass_worker_limit") <- config$worker_limit
  param
}

.rc_parallel_param_for_tasks <- function(BPPARAM = NULL, n_tasks) {
  n_tasks <- suppressWarnings(as.integer(n_tasks[[1L]]))
  if (!is.finite(n_tasks) || n_tasks <= 1L || identical(BPPARAM, FALSE)) {
    return(FALSE)
  }
  if (is.null(BPPARAM)) {
    BPPARAM <- rc_default_bpparam()
    if (is.null(BPPARAM)) return(FALSE)
  }
  if (!requireNamespace("BiocParallel", quietly = TRUE) ||
      !methods::is(BPPARAM, "BiocParallelParam")) {
    stop("`BPPARAM` must be NULL, FALSE, or a BiocParallelParam object.",
         call. = FALSE)
  }
  limit <- .rc_bpparam_worker_limit(BPPARAM)
  effective <- min(limit, n_tasks)
  if (effective <= 1L) return(FALSE)

  if (isTRUE(BiocParallel::bpisup(BPPARAM))) {
    return(BPPARAM)
  }
  current <- max(1L, as.integer(BiocParallel::bpnworkers(BPPARAM)))
  if (current == effective) {
    attr(BPPARAM, "regcompass_worker_limit") <- limit
    attr(BPPARAM, "regcompass_effective_workers") <- effective
    return(BPPARAM)
  }
  backend <- .rc_bpparam_backend(BPPARAM)
  tuned <- rc_default_bpparam(workers = effective, backend = backend)
  if (is.null(tuned)) return(FALSE)
  attr(tuned, "regcompass_worker_limit") <- limit
  attr(tuned, "regcompass_effective_workers") <- effective
  tuned
}

#' Apply a function with optional BiocParallel control
#'
#' Every task, including tasks submitted to a caller-started pool, establishes a
#' one-thread numerical/solver environment inside the worker. Package-managed
#' pools are sized as `min(number of independent tasks, worker cap)`, stopped
#' after the dispatch, and followed by full garbage collection.
#'
#' @param X A vector or list.
#' @param FUN Function applied to each element.
#' @param BPPARAM Internal `NULL`, `FALSE`, or a `BiocParallelParam` template.
#' @param ... Additional arguments.
#' @return A list.
rc_parallel_lapply <- function(X, FUN, BPPARAM = NULL, ...) {
  if (!is.function(FUN)) stop("`FUN` must be a function.", call. = FALSE)
  extra <- list(...)
  caller_libpaths <- .libPaths()
  worker_fun <- function(x) {
    .libPaths(unique(c(caller_libpaths, .libPaths())))
    .rc_with_internal_single_thread(function() {
      do.call(FUN, c(list(x), extra))
    })
  }
  if (identical(BPPARAM, FALSE) || length(X) <= 1L) {
    return(lapply(X, worker_fun))
  }
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

  original <- BPPARAM
  BPPARAM <- .rc_parallel_param_for_tasks(BPPARAM, length(X))
  if (identical(BPPARAM, FALSE) || is.null(BPPARAM)) {
    return(lapply(X, worker_fun))
  }
  was_started <- isTRUE(BiocParallel::bpisup(BPPARAM))
  thread_state <- NULL
  if (!was_started) {
    thread_state <- .rc_set_internal_single_thread()
    on.exit({
      .rc_release_bpparam(BPPARAM)
      .rc_restore_internal_threads(thread_state)
      invisible(gc(verbose = FALSE, full = TRUE))
    }, add = TRUE)
    BiocParallel::bpstart(BPPARAM)
  } else if (!identical(original, BPPARAM)) {
    stop("Internal error: a started worker pool cannot be resized.",
         call. = FALSE)
  }
  BiocParallel::bplapply(X, worker_fun, BPPARAM = BPPARAM)
}
