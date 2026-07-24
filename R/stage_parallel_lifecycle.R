.rc_internal_thread_env <- function() {
  c(
    OMP_NUM_THREADS = "1",
    OMP_DYNAMIC = "FALSE",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    MKL_DYNAMIC = "FALSE",
    VECLIB_MAXIMUM_THREADS = "1",
    BLIS_NUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1",
    RCPP_PARALLEL_NUM_THREADS = "1"
  )
}

.rc_set_internal_single_thread <- function() {
  desired <- .rc_internal_thread_env()
  old_env <- Sys.getenv(names(desired), unset = NA_character_)
  old_options <- options(
    mc.cores = 1L,
    RegCompassR.internal_workers = 1L
  )
  do.call(Sys.setenv, as.list(desired))
  list(env = old_env, options = old_options)
}

.rc_restore_internal_threads <- function(state) {
  if (!is.list(state)) return(invisible(NULL))
  old_env <- state$env %||% character()
  if (length(old_env)) {
    missing <- names(old_env)[is.na(old_env)]
    present <- names(old_env)[!is.na(old_env)]
    if (length(present)) {
      do.call(Sys.setenv, as.list(old_env[present]))
    }
    if (length(missing)) Sys.unsetenv(missing)
  }
  if (is.list(state$options)) do.call(options, state$options)
  invisible(NULL)
}

.rc_with_internal_single_thread <- function(FUN) {
  if (!is.function(FUN)) stop("`FUN` must be a function.", call. = FALSE)
  state <- .rc_set_internal_single_thread()
  on.exit(.rc_restore_internal_threads(state), add = TRUE)
  FUN()
}

.rc_stage_worker_config <- function(workers, argument = "workers") {
  if (length(workers) != 1L || is.na(workers) || !is.finite(workers)) {
    stop("`", argument, "` must be one positive integer.", call. = FALSE)
  }
  workers <- suppressWarnings(as.integer(workers))
  if (workers < 1L) {
    stop("`", argument, "` must be at least 1.", call. = FALSE)
  }
  rc_parallel_config(workers = workers, backend = "auto")
}

.rc_with_stage_workers <- function(workers, FUN, argument = "workers") {
  if (!is.function(FUN)) stop("`FUN` must be a function.", call. = FALSE)
  config <- .rc_stage_worker_config(workers, argument = argument)
  thread_state <- .rc_set_internal_single_thread()
  on.exit(.rc_restore_internal_threads(thread_state), add = TRUE)

  param <- .rc_phase_bpparam(config$workers, backend = "auto")
  on.exit({
    if (!identical(param, FALSE) && !is.null(param)) {
      .rc_release_bpparam(param)
    }
    param <- FALSE
    invisible(gc(verbose = FALSE, full = TRUE))
  }, add = TRUE)

  if (!identical(param, FALSE) && !is.null(param)) {
    BiocParallel::bpstart(param)
  }

  FUN(param, config)
}
