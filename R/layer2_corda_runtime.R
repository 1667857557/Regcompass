# Runtime context and parallel scheduling for exact pinned Python CORDA2.

.rc_layer2_completion_context <- new.env(parent = emptyenv())
.rc_layer2_completion_context$active <- FALSE
.rc_layer2_completion_context$model_completion <- "fastcore"
.rc_layer2_completion_context$reaction_evidence <- NULL
.rc_layer2_completion_context$corda_options <- NULL

.rc_is_corda2_options <- function(options) {
  is.list(options) && identical(
    as.character(options$algorithm %||% ""),
    "resendislab_python_CORDA2_c02e06d_exact_semantics"
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
  setter <- get0(
    "bptasks<-", envir = asNamespace("BiocParallel"),
    mode = "function", inherits = FALSE
  )
  if (is.function(setter)) {
    tuned <- tryCatch(
      setter(BPPARAM, n_tasks),
      error = function(e) BPPARAM
    )
  } else {
    tuned <- BPPARAM
  }
  attr(tuned, "regcompass_corda2_dynamic_tasks") <- n_tasks
  tuned
}

.rc_corda_should_outer_parallel <- function(n_tasks, pool_workers) {
  as.integer(n_tasks) > 1L && as.integer(pool_workers) > 1L
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
