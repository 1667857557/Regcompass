# Layer 2 parallel context helpers.

.rc_layer2_parallel_context <- new.env(parent = emptyenv())
.rc_layer2_parallel_context$active <- FALSE
.rc_layer2_parallel_context$parallel <- TRUE
.rc_layer2_parallel_context$BPPARAM <- NULL
.rc_layer2_parallel_context$nested_serial <- FALSE

.rc_layer2_enter_parallel_context <- function(parallel, BPPARAM) {
  previous <- as.list(.rc_layer2_parallel_context)
  .rc_layer2_parallel_context$active <- TRUE
  .rc_layer2_parallel_context$parallel <- isTRUE(parallel)
  .rc_layer2_parallel_context$BPPARAM <- BPPARAM
  .rc_layer2_parallel_context$nested_serial <- FALSE
  previous
}

.rc_layer2_restore_parallel_context <- function(previous) {
  rm(
    list = ls(.rc_layer2_parallel_context, all.names = TRUE),
    envir = .rc_layer2_parallel_context
  )
  list2env(previous, envir = .rc_layer2_parallel_context)
  invisible(NULL)
}

.rc_atomic_save_rds <- function(object, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(
    pattern = paste0(basename(file), ".tmp_"),
    tmpdir = dirname(file)
  )
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  saveRDS(object, temporary)
  if (!file.rename(temporary, file)) {
    if (!file.copy(temporary, file, overwrite = TRUE)) {
      stop("Cannot persist Layer 2 task result: ", file, call. = FALSE)
    }
    unlink(temporary, force = TRUE)
  }
  invisible(file)
}

.rc_layer2_task_bpparam <- function() {
  if (!isTRUE(.rc_layer2_parallel_context$active) ||
      !isTRUE(.rc_layer2_parallel_context$parallel) ||
      isTRUE(.rc_layer2_parallel_context$nested_serial)) {
    return(FALSE)
  }
  .rc_layer2_parallel_context$BPPARAM
}

.rc_layer2_pool_workers <- function(BPPARAM) {
  if (identical(BPPARAM, FALSE) || is.null(BPPARAM)) return(1L)
  recorded <- suppressWarnings(as.integer(
    attr(BPPARAM, "regcompass_layer2_effective_workers") %||% NA_integer_
  ))
  if (is.finite(recorded) && recorded >= 1L) return(recorded)
  if (!is.null(BPPARAM) &&
      requireNamespace("BiocParallel", quietly = TRUE) &&
      methods::is(BPPARAM, "BiocParallelParam")) {
    return(max(1L, as.integer(BiocParallel::bpnworkers(BPPARAM))))
  }
  config <- rc_parallel_config(workers = NULL, backend = "auto")
  max(1L, as.integer(config$workers))
}

.rc_layer2_tune_task_bpparam <- function(BPPARAM, n_tasks) {
  n_tasks <- max(1L, as.integer(n_tasks[[1L]]))
  if (identical(BPPARAM, FALSE) || n_tasks <= 1L) return(BPPARAM)

  if (is.null(BPPARAM)) {
    requested <- max(1L, as.integer(rc_available_workers(default = 1L)))
    BPPARAM <- rc_default_bpparam(workers = min(requested, n_tasks))
    if (is.null(BPPARAM)) return(NULL)
  } else {
    requested <- .rc_layer2_pool_workers(BPPARAM)
  }

  if (!requireNamespace("BiocParallel", quietly = TRUE) ||
      !methods::is(BPPARAM, "BiocParallelParam")) {
    attr(BPPARAM, "regcompass_layer2_requested_workers") <- requested
    attr(BPPARAM, "regcompass_layer2_effective_workers") <-
      min(requested, n_tasks)
    return(BPPARAM)
  }

  tuned <- BPPARAM
  started <- isTRUE(BiocParallel::bpisup(tuned))
  effective <- min(
    n_tasks,
    max(1L, as.integer(BiocParallel::bpnworkers(tuned)))
  )
  if (!started && effective < BiocParallel::bpnworkers(tuned)) {
    worker_setter <- get0(
      "bpnworkers<-", envir = asNamespace("BiocParallel"),
      mode = "function", inherits = FALSE
    )
    if (is.function(worker_setter)) {
      tuned <- tryCatch(
        worker_setter(tuned, effective),
        error = function(e) tuned
      )
    }
  }

  task_setter <- get0(
    "bptasks<-", envir = asNamespace("BiocParallel"),
    mode = "function", inherits = FALSE
  )
  if (is.function(task_setter)) {
    tuned <- tryCatch(
      task_setter(tuned, n_tasks),
      error = function(e) tuned
    )
  }

  effective <- min(
    n_tasks,
    max(1L, as.integer(BiocParallel::bpnworkers(tuned)))
  )
  attr(tuned, "regcompass_layer2_requested_workers") <- requested
  attr(tuned, "regcompass_layer2_effective_workers") <- effective
  attr(tuned, "regcompass_layer2_dynamic_tasks") <- n_tasks
  tuned
}

.rc_layer2_should_outer_parallel <- function(n_tasks, pool_workers) {
  as.integer(n_tasks) > 1L && as.integer(pool_workers) > 1L
}
