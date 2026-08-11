# Final budget-aware dispatcher.
#
# The outer BiocParallel pool is resized to the number of immediately runnable
# jobs, but nested Pando budgets are computed from the original RegCompass worker
# cap rather than from that resized outer pool. For N cell-type jobs and cap W,
# the per-GRN budgets sum to W when N <= W. Each GRN budget includes its outer
# worker, so nested target pools may use at most budget - 1 processes.

.rc_is_hierarchical_pando_dispatch <- function(FUN, budgets, task_keys) {
  if (!length(task_keys) || length(budgets) <= 1L || !any(budgets > 1L)) {
    return(FALSE)
  }
  condition_fun <- get0(
    ".rc_condition_multitask_fit_task", envir = asNamespace("RegCompassR"),
    inherits = FALSE, ifnotfound = NULL
  )
  standard_fun <- get0(
    ".rc_run_standard_pando_celltype_job", envir = asNamespace("RegCompassR"),
    inherits = FALSE, ifnotfound = NULL
  )
  identical(FUN, condition_fun) || identical(FUN, standard_fun)
}

.rc_hierarchical_outer_param <- function(BPPARAM, workers, global_limit) {
  workers <- max(1L, as.integer(workers[[1L]]))
  if (workers <= 1L) return(FALSE)
  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    stop("BiocParallel is required for hierarchical Pando execution.",
         call. = FALSE)
  }
  param <- BiocParallel::SnowParam(
    workers = workers,
    type = "SOCK",
    progressbar = FALSE,
    exportglobals = TRUE,
    exportvariables = TRUE
  )
  attr(param, "regcompass_worker_limit") <- as.integer(global_limit)
  attr(param, "regcompass_effective_workers") <- workers
  attr(param, "regcompass_hierarchical_outer_pool") <- TRUE
  param
}

rc_parallel_lapply <- function(X, FUN, BPPARAM = NULL, ...) {
  if (!length(X)) return(list())
  if (!is.function(FUN)) stop("`FUN` must be a function.", call. = FALSE)
  extra <- list(...)
  caller_libpaths <- .libPaths()
  internal_env <- .rc_internal_thread_env()

  validate_param <- function(param) {
    if (is.null(param) || identical(param, FALSE)) return(invisible(TRUE))
    if (is.logical(param) || !requireNamespace("BiocParallel", quietly = TRUE) ||
        !methods::is(param, "BiocParallelParam")) {
      stop("`BPPARAM` must be NULL, FALSE, or a BiocParallelParam object.",
           call. = FALSE)
    }
    invisible(TRUE)
  }
  validate_param(BPPARAM)

  global_limit <- if (!is.null(BPPARAM) && !identical(BPPARAM, FALSE)) {
    .rc_bpparam_worker_limit(BPPARAM, default = 1L)
  } else {
    1L
  }
  global_backend <- if (!is.null(BPPARAM) && !identical(BPPARAM, FALSE)) {
    .rc_bpparam_backend(BPPARAM)
  } else {
    "serial"
  }
  budgets <- .rc_allocate_condition_target_workers(global_limit, length(X))
  task_keys <- if (all(vapply(X, function(x) {
    is.list(x) && is.character(x$cell_type) && length(x$cell_type) == 1L &&
      !is.na(x$cell_type) && nzchar(x$cell_type)
  }, logical(1)))) {
    vapply(X, `[[`, character(1), "cell_type")
  } else {
    character()
  }
  if (length(task_keys) == length(X) && !anyDuplicated(task_keys)) {
    names(budgets) <- task_keys
  }
  hierarchical_pando <- .rc_is_hierarchical_pando_dispatch(
    FUN, budgets, task_keys
  )

  worker_fun <- function(x) {
    .libPaths(unique(c(caller_libpaths, .libPaths())))
    task_workers <- 1L
    if (length(task_keys) == length(X) && !anyDuplicated(task_keys) &&
        is.list(x) && is.character(x$cell_type) &&
        length(x$cell_type) == 1L && x$cell_type %in% task_keys) {
      task_workers <- as.integer(budgets[[x$cell_type]])
    }
    dispatch_env <- c(
      REGCOMPASS_TASK_WORKERS = as.character(task_workers),
      REGCOMPASS_DISPATCH_WORKER_LIMIT = as.character(global_limit),
      REGCOMPASS_DISPATCH_TASKS = as.character(length(X)),
      REGCOMPASS_DISPATCH_BACKEND = as.character(
        if (hierarchical_pando) "snow_hierarchical" else global_backend
      )
    )
    env_all <- c(internal_env, dispatch_env)
    old_env <- Sys.getenv(names(env_all), unset = NA_character_)
    old_options <- base::options(
      mc.cores = 1L,
      RegCompassR.internal_workers = 1L
    )
    do.call(Sys.setenv, as.list(env_all))
    on.exit({
      missing <- names(old_env)[is.na(old_env)]
      present <- names(old_env)[!is.na(old_env)]
      if (length(present)) do.call(Sys.setenv, as.list(old_env[present]))
      if (length(missing)) Sys.unsetenv(missing)
      do.call(base::options, old_options)
      invisible(gc(verbose = FALSE, full = TRUE))
    }, add = TRUE)
    do.call(FUN, c(list(x), extra))
  }

  if (identical(BPPARAM, FALSE) || is.null(BPPARAM)) {
    return(lapply(X, worker_fun))
  }
  if (length(X) == 1L) {
    return(lapply(X, worker_fun))
  }

  outer_workers <- min(global_limit, length(X))
  outer_param <- if (hierarchical_pando) {
    .rc_hierarchical_outer_param(BPPARAM, outer_workers, global_limit)
  } else {
    .rc_parallel_param_for_tasks(BPPARAM, length(X))
  }
  if (identical(outer_param, FALSE) || is.null(outer_param)) {
    return(lapply(X, worker_fun))
  }
  was_started <- isTRUE(BiocParallel::bpisup(outer_param))
  thread_state <- NULL
  if (!was_started) {
    thread_state <- .rc_set_internal_single_thread()
    BiocParallel::bpstart(outer_param)
  }
  on.exit({
    if (!was_started) {
      try(BiocParallel::bpstop(outer_param), silent = TRUE)
      .rc_restore_internal_threads(thread_state)
    }
    invisible(gc(verbose = FALSE, full = TRUE))
  }, add = TRUE)
  BiocParallel::bplapply(X, worker_fun, BPPARAM = outer_param)
}
