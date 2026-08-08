#' Detect the RegCompass worker budget
#'
#' Explicit RegCompass settings take precedence. Without an explicit setting,
#' RegCompass uses at most `default` workers and never intentionally exceeds the
#' CPU allocation detected from the scheduler or current machine.
#'
#' @param default Fallback worker upper bound. The package default is 10.
#' @return A positive integer worker count.
rc_available_workers <- function(default = 10L) {
  explicit <- c(
    getOption("RegCompassR.workers", NA),
    Sys.getenv("REGCOMPASS_WORKERS", unset = NA_character_),
    Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA_character_),
    Sys.getenv("NSLOTS", unset = NA_character_)
  )
  explicit <- suppressWarnings(as.integer(explicit))
  explicit <- explicit[is.finite(explicit) & explicit >= 1L]
  if (length(explicit)) return(max(1L, explicit[[1L]]))

  default <- suppressWarnings(as.integer(default[[1L]]))
  if (!is.finite(default) || default < 1L) default <- 10L
  available <- NA_integer_
  if (requireNamespace("future", quietly = TRUE)) {
    available <- tryCatch(
      suppressWarnings(as.integer(future::availableCores()[[1L]])),
      error = function(e) NA_integer_
    )
  }
  if (!is.finite(available) || available < 1L) {
    cores <- suppressWarnings(as.integer(parallel::detectCores(logical = TRUE)[[1L]]))
    if (is.finite(cores) && cores >= 1L) available <- max(1L, cores - 1L)
  }
  if (!is.finite(available) || available < 1L) return(default)
  max(1L, min(default, available))
}

.rc_normalize_worker_budget <- function(workers = NULL, argument = "workers") {
  if (is.null(workers)) workers <- rc_available_workers(default = 10L)
  if (length(workers) != 1L || is.na(workers) || !is.finite(workers)) {
    stop("`", argument, "` must be one positive integer or NULL.", call. = FALSE)
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

#' Resolve the platform-aware RegCompass parallel budget
#'
#' `workers` is the only user-facing parallel budget. The default budget is 10.
#' `backend = "auto"` selects a SOCK cluster on Windows and forked multicore
#' workers on Linux/macOS. Individual operations automatically use no more than
#' `min(number_of_independent_tasks, workers)` workers. Sequential execution is
#' used when one worker is requested or BiocParallel is unavailable.
#'
#' @param workers Global worker upper bound. `NULL` uses an explicit
#' `options(RegCompassR.workers)`/environment setting when present, otherwise a
#' machine-aware maximum of 10 workers.
#' @param backend Requested backend.
#' @return A list describing requested and resolved execution settings.
#' @export
rc_parallel_config <- function(
    workers = NULL,
    backend = c("auto", "serial", "snow", "multicore")) {
  backend <- match.arg(backend)
  requested_workers <- workers
  workers <- .rc_normalize_worker_budget(workers)
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
    worker_budget = workers,
    biocparallel_available = available
  )
}

#' Build a platform-aware RegCompass parallel backend
#'
#' This is an internal-facing backend constructor. Workflow users normally set
#' only a `workers` budget on the public workflow or step function. The backend's
#' task-level progress bar follows `options(RegCompassR.progress = TRUE/FALSE)`.
#'
#' @param workers Worker upper bound for this backend.
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
  attr(param, "regcompass_worker_budget") <- config$worker_budget
  param
}

.rc_tune_bpparam_to_tasks <- function(BPPARAM, n_tasks, workers = NULL) {
  n_tasks <- max(1L, as.integer(n_tasks[[1L]]))
  budget <- .rc_normalize_worker_budget(
    workers %||% attr(BPPARAM, "regcompass_worker_budget") %||%
      rc_available_workers(default = 10L)
  )
  if (identical(BPPARAM, FALSE)) return(FALSE)
  if (is.null(BPPARAM)) {
    return(.rc_task_bpparam(workers = budget, n_tasks = n_tasks))
  }
  if (!requireNamespace("BiocParallel", quietly = TRUE) ||
      !methods::is(BPPARAM, "BiocParallelParam")) return(BPPARAM)

  effective <- min(
    n_tasks, budget,
    max(1L, as.integer(BiocParallel::bpnworkers(BPPARAM)))
  )
  tuned <- BPPARAM
  if (!isTRUE(BiocParallel::bpisup(tuned)) &&
      effective < BiocParallel::bpnworkers(tuned)) {
    setter <- get0(
      "bpnworkers<-", envir = asNamespace("BiocParallel"),
      mode = "function", inherits = FALSE
    )
    if (is.function(setter)) {
      tuned <- tryCatch(setter(tuned, effective), error = function(e) tuned)
    }
  }
  task_setter <- get0(
    "bptasks<-", envir = asNamespace("BiocParallel"),
    mode = "function", inherits = FALSE
  )
  if (is.function(task_setter)) {
    tuned <- tryCatch(task_setter(tuned, n_tasks), error = function(e) tuned)
  }
  attr(tuned, "regcompass_worker_budget") <- budget
  attr(tuned, "regcompass_effective_workers") <- effective
  attr(tuned, "regcompass_dynamic_tasks") <- n_tasks
  tuned
}

.rc_task_bpparam <- function(
    workers = NULL, n_tasks = NULL,
    backend = c("auto", "serial", "snow", "multicore")) {
  backend <- match.arg(backend)
  budget <- .rc_normalize_worker_budget(workers)
  if (is.null(n_tasks)) n_tasks <- budget
  n_tasks <- max(1L, as.integer(n_tasks[[1L]]))
  effective <- min(budget, n_tasks)
  param <- rc_default_bpparam(workers = effective, backend = backend)
  if (is.null(param)) return(FALSE)
  attr(param, "regcompass_worker_budget") <- budget
  attr(param, "regcompass_effective_workers") <- effective
  attr(param, "regcompass_dynamic_tasks") <- n_tasks
  param
}

#' Apply a function under the RegCompass worker budget
#'
#' Each call automatically reduces the backend to at most the number of
#' independent tasks. Every task establishes a one-thread numerical/solver
#' environment inside the worker, preventing nested BLAS/HiGHS oversubscription.
#' Package-managed pools are stopped and followed by full garbage collection.
#'
#' @param X A vector or list.
#' @param FUN Function applied to each element.
#' @param BPPARAM Internal `BiocParallelParam`, `NULL`, or `FALSE`.
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
  BPPARAM <- .rc_tune_bpparam_to_tasks(BPPARAM, length(X))
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
  }
  BiocParallel::bplapply(X, worker_fun, BPPARAM = BPPARAM)
}
