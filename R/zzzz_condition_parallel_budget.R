# Worker-budget helpers for bounded hierarchical Pando scheduling.

.rc_allocate_condition_target_workers <- function(worker_limit, n_tasks) {
  worker_limit <- suppressWarnings(as.integer(worker_limit[[1L]]))
  n_tasks <- suppressWarnings(as.integer(n_tasks[[1L]]))
  if (!is.finite(worker_limit) || worker_limit < 1L ||
      !is.finite(n_tasks) || n_tasks < 1L) {
    stop("Worker limit and task count must be positive integers.", call. = FALSE)
  }
  if (n_tasks >= worker_limit) return(rep.int(1L, n_tasks))
  base <- worker_limit %/% n_tasks
  remainder <- worker_limit %% n_tasks
  as.integer(base + (seq_len(n_tasks) <= remainder))
}

.rc_condition_nested_target_bpparam <- function(workers) {
  workers <- suppressWarnings(as.integer(workers[[1L]]))
  if (!is.finite(workers) || workers <= 1L) return(NULL)
  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    stop("BiocParallel is required for Pando target parallelism.",
         call. = FALSE)
  }
  # Always use an isolated SOCK pool for nested target work. Forking from an
  # already parallel cell-type worker is deliberately avoided on Unix too.
  param <- BiocParallel::SnowParam(
    workers = workers,
    type = "SOCK",
    progressbar = FALSE,
    exportglobals = TRUE,
    exportvariables = TRUE
  )
  attr(param, "regcompass_worker_limit") <- workers
  attr(param, "regcompass_effective_workers") <- workers
  param
}
