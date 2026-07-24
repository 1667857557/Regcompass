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
  param <- .rc_phase_bpparam(config$workers, backend = "auto")
  started <- FALSE

  if (!identical(param, FALSE) && !is.null(param)) {
    BiocParallel::bpstart(param)
    started <- TRUE
  }

  on.exit({
    if (started || (!identical(param, FALSE) && !is.null(param))) {
      .rc_release_bpparam(param)
    }
    param <- FALSE
    invisible(gc(verbose = FALSE, full = TRUE))
  }, add = TRUE)

  FUN(param, config)
}
