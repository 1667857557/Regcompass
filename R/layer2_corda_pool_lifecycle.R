# Keep one worker pool alive across independent CORDA2 model instances.

.rc_layer2_requested_corda2 <- function(layer2_args) {
  model_params <- if (is.list(layer2_args)) {
    layer2_args$model_params %||% list()
  } else {
    list()
  }
  requested <- as.character(model_params$model_completion %||% "fastcore")
  length(requested) == 1L && !is.na(requested) &&
    requested %in% c("corda2", "corda")
}

.rc_prepare_corda_worker_pool <- function(
    layer2_args, parallel = TRUE, BPPARAM = NULL) {
  state <- list(
    is_corda2 = .rc_layer2_requested_corda2(layer2_args),
    BPPARAM = BPPARAM,
    origin = "not_used",
    started_here = FALSE,
    thread_state = NULL
  )
  if (!isTRUE(state$is_corda2) || !isTRUE(parallel)) return(state)
  if (is.null(state$BPPARAM)) {
    state$BPPARAM <- rc_default_bpparam()
    state$origin <- if (is.null(state$BPPARAM)) {
      "serial_fallback"
    } else {
      "package_default"
    }
  } else {
    state$origin <- "caller_supplied"
  }
  if (is.null(state$BPPARAM)) return(state)
  if (!requireNamespace("BiocParallel", quietly = TRUE) ||
      !methods::is(state$BPPARAM, "BiocParallelParam")) {
    stop(
      "CORDA2 outer-model parallel execution requires a ",
      "BiocParallelParam object.",
      call. = FALSE
    )
  }
  if (!isTRUE(BiocParallel::bpisup(state$BPPARAM))) {
    state$thread_state <- .rc_set_internal_single_thread()
    start_error <- tryCatch({
      BiocParallel::bpstart(state$BPPARAM)
      NULL
    }, error = function(e) e)
    if (inherits(start_error, "error")) {
      .rc_restore_internal_threads(state$thread_state)
      stop(
        "Unable to start the CORDA2 BiocParallel worker pool: ",
        conditionMessage(start_error),
        call. = FALSE
      )
    }
    state$started_here <- TRUE
  }
  state
}

.rc_release_corda_worker_pool <- function(state) {
  if (!is.list(state)) return(invisible(NULL))
  if (isTRUE(state$started_here) &&
      requireNamespace("BiocParallel", quietly = TRUE) &&
      methods::is(state$BPPARAM, "BiocParallelParam") &&
      isTRUE(BiocParallel::bpisup(state$BPPARAM))) {
    .rc_release_bpparam(state$BPPARAM)
  }
  if (!is.null(state$thread_state)) {
    .rc_restore_internal_threads(state$thread_state)
  }
  invisible(gc(verbose = FALSE, full = TRUE))
  invisible(NULL)
}
