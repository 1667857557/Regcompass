# Runtime context and parallel scheduling for original MATLAB CORDA2.

.rc_layer2_completion_context <- new.env(parent = emptyenv())
.rc_layer2_completion_context$active <- FALSE
.rc_layer2_completion_context$model_completion <- "corda2"
.rc_layer2_completion_context$reaction_evidence <- NULL
.rc_layer2_completion_context$corda_options <- NULL

.rc_is_corda2_options <- function(options) {
  is.list(options) && identical(
    as.character(options$algorithm %||% ""),
    "schultzdre_MATLAB_CORDA2_original_semantics"
  )
}

.rc_corda_tune_task_bpparam <- function(BPPARAM, n_tasks) {
  n_tasks <- max(1L, as.integer(n_tasks[[1L]]))
  if (identical(BPPARAM, FALSE) || is.null(BPPARAM) || n_tasks <= 1L) {
    return(BPPARAM)
  }
  if (!requireNamespace("BiocParallel", quietly = TRUE) ||
      !methods::is(BPPARAM, "BiocParallelParam")) {
    return(BPPARAM)
  }

  # BiocParallelParam is a reference class. Never write `bptasks` into the
  # shared Layer 2 template: a small structural task count would otherwise
  # persist into later Vmax and reaction-level Step 2 scoring dispatches.
  tuned <- .rc_parallel_param_for_tasks(BPPARAM, n_tasks)
  if (identical(tuned, FALSE) || is.null(tuned)) return(tuned)
  setter <- get0(
    "bptasks<-", envir = asNamespace("BiocParallel"),
    mode = "function", inherits = FALSE
  )
  if (is.function(setter)) {
    tuned <- tryCatch(
      setter(tuned, n_tasks),
      error = function(e) tuned
    )
  }
  attr(tuned, "regcompass_corda2_dynamic_tasks") <- n_tasks
  tuned
}

.rc_corda_stage_parallel_requested <- function() {
  context <- get0(
    ".rc_layer2_parallel_context", mode = "environment", inherits = TRUE
  )
  if (!is.environment(context) || !isTRUE(context$active) ||
      !isTRUE(context$parallel)) return(FALSE)
  template <- context$BPPARAM
  if (identical(template, FALSE)) return(FALSE)
  if (!is.null(template) &&
      requireNamespace("BiocParallel", quietly = TRUE) &&
      methods::is(template, "BiocParallelParam")) {
    return(as.integer(BiocParallel::bpnworkers(template)) > 1L)
  }
  available <- rc_available_workers(default = 1L)
  requireNamespace("BiocParallel", quietly = TRUE) &&
    as.integer(available) > 1L
}

.rc_corda_should_outer_parallel <- function(n_tasks, pool_workers) {
  FALSE
}

.rc_corda_pool_workers <- function(BPPARAM) {
  if (identical(BPPARAM, FALSE)) return(1L)
  if (!is.null(BPPARAM) &&
      requireNamespace("BiocParallel", quietly = TRUE) &&
      methods::is(BPPARAM, "BiocParallelParam")) {
    return(max(1L, BiocParallel::bpnworkers(BPPARAM)))
  }
  config <- rc_parallel_config(workers = NULL, backend = "auto")
  max(1L, config$workers)
}
