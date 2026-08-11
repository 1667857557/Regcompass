# Bounded hierarchical scheduling for condition-GRN cell types and Pando targets.

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

.rc_parallel_lapply <- function(X, FUN, BPPARAM = NULL, ...) {
  if (!is.function(FUN)) stop("`FUN` must be a function.", call. = FALSE)
  extra <- list(...)
  caller_libpaths <- .libPaths()
  internal_env <- .rc_internal_thread_env()
  task_keys <- if (length(X) && all(vapply(X, function(x) {
    is.list(x) && is.character(x$cell_type) && length(x$cell_type) == 1L &&
      !is.na(x$cell_type) && nzchar(x$cell_type)
  }, logical(1)))) {
    vapply(X, `[[`, character(1), "cell_type")
  } else {
    character()
  }

  make_worker_fun <- function(dispatch_limit, dispatch_backend) {
    budgets <- .rc_allocate_condition_target_workers(dispatch_limit, length(X))
    if (length(task_keys) == length(X) && !anyDuplicated(task_keys)) {
      names(budgets) <- task_keys
    }
    function(x) {
      .libPaths(unique(c(caller_libpaths, .libPaths())))
      task_workers <- 1L
      if (length(task_keys) == length(X) && !anyDuplicated(task_keys) &&
          is.list(x) && is.character(x$cell_type) &&
          length(x$cell_type) == 1L && x$cell_type %in% task_keys) {
        task_workers <- as.integer(budgets[[x$cell_type]])
      }
      dispatch_env <- c(
        REGCOMPASS_TASK_WORKERS = as.character(task_workers),
        REGCOMPASS_DISPATCH_WORKER_LIMIT = as.character(dispatch_limit),
        REGCOMPASS_DISPATCH_TASKS = as.character(length(X)),
        REGCOMPASS_DISPATCH_BACKEND = as.character(dispatch_backend)
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
      }, add = TRUE)
      do.call(FUN, c(list(x), extra))
    }
  }

  if (identical(BPPARAM, FALSE) || length(X) <= 1L) {
    worker_fun <- make_worker_fun(1L, "serial")
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
      stop("`BPPARAM` must be NULL, FALSE, or a BiocParallelParam object.",
           call. = FALSE)
    }
  }

  original <- BPPARAM
  BPPARAM <- .rc_parallel_param_for_tasks(BPPARAM, length(X))
  if (identical(BPPARAM, FALSE) || is.null(BPPARAM)) {
    worker_fun <- make_worker_fun(1L, "serial")
    return(lapply(X, worker_fun))
  }
  dispatch_limit <- .rc_bpparam_worker_limit(BPPARAM, default = 1L)
  dispatch_backend <- .rc_bpparam_backend(BPPARAM)
  worker_fun <- make_worker_fun(dispatch_limit, dispatch_backend)
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

.rc_condition_nested_target_bpparam <- function(workers) {
  workers <- suppressWarnings(as.integer(workers[[1L]]))
  if (!is.finite(workers) || workers <= 1L) return(NULL)
  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    stop("BiocParallel is required for condition target parallelism.",
         call. = FALSE)
  }
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

.rc_condition_multitask_fit_task <- function(
    task, target_genes, condition_col, celltype_col, min_cells,
    pando_infer_args, inner_parallel = FALSE, PANDO_BPPARAM = NULL) {
  if (!is.list(task) || !inherits(task$grn, "GRNData") ||
      !is.character(task$cell_type) || length(task$cell_type) != 1L) {
    stop("Invalid condition-GRN multi-task fit task.", call. = FALSE)
  }
  ridge_control <- pando_infer_args$condition_ridge_control %||% list()
  threshold <- suppressWarnings(as.numeric(pando_infer_args$padj_threshold))
  if (length(threshold) != 1L || !is.finite(threshold) ||
      threshold <= 0 || threshold >= 1) {
    stop("Condition Pando padj_threshold must be in (0, 1).",
         call. = FALSE)
  }

  allocated_workers <- suppressWarnings(as.integer(
    Sys.getenv("REGCOMPASS_TASK_WORKERS", unset = "1")
  ))
  if (!is.finite(allocated_workers) || allocated_workers < 1L) {
    allocated_workers <- 1L
  }
  target_param <- NULL
  target_backend <- "serial"
  if (isTRUE(inner_parallel) && !is.null(PANDO_BPPARAM) &&
      !identical(PANDO_BPPARAM, FALSE)) {
    target_param <- PANDO_BPPARAM
    allocated_workers <- .rc_bpparam_worker_limit(PANDO_BPPARAM, default = 1L)
    target_backend <- .rc_bpparam_backend(PANDO_BPPARAM)
  } else if (allocated_workers > 1L) {
    target_param <- .rc_condition_nested_target_bpparam(allocated_workers)
    target_backend <- "snow"
  }
  target_parallel <- allocated_workers > 1L && !is.null(target_param)

  old_r_libs <- Sys.getenv("R_LIBS", unset = NA_character_)
  if (target_parallel) {
    Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
    on.exit({
      if (is.na(old_r_libs)) Sys.unsetenv("R_LIBS") else Sys.setenv(R_LIBS = old_r_libs)
    }, add = TRUE)
  }

  args <- list(
    object = task$grn,
    cell_type_col = celltype_col,
    condition_col = condition_col,
    cell_type = task$cell_type,
    genes = target_genes,
    network_name = "regcompass_condition_grn",
    rna_layer = pando_infer_args$rna_layer %||% "data",
    peak_layer = pando_infer_args$peak_layer %||% "data",
    peak_value_type = pando_infer_args$peak_value_type %||% "normalized",
    tf_cor = pando_infer_args$tf_cor,
    peak_cor = pando_infer_args$peak_cor,
    min_cells_per_condition = as.integer(min_cells),
    small_condition_action = "error",
    adjust_method = "BH",
    padj_threshold = threshold,
    rank_action = pando_infer_args$rank_action,
    min_residual_df = pando_infer_args$min_residual_df,
    parallel = target_parallel,
    parallel_scope = "target",
    overwrite = TRUE,
    fallback_args = list(condition_ridge_control = ridge_control),
    verbose = FALSE
  )
  if (target_parallel) args$BPPARAM <- target_param
  fitted <- do.call(Pando::infer_condition_grn, args)
  fits <- Pando::condition_grn_fit(fitted)
  if (inherits(fits, "ConditionGRNFit")) {
    fit <- fits
  } else if (is.list(fits) && task$cell_type %in% names(fits)) {
    fit <- fits[[task$cell_type]]
  } else if (is.list(fits) && length(fits) == 1L &&
             inherits(fits[[1L]], "ConditionGRNFit")) {
    fit <- fits[[1L]]
  } else {
    stop("Pando multi-task condition fit was not returned for cell type `",
         task$cell_type, "`.", call. = FALSE)
  }
  if (!identical(as.character(fit$cell_type), task$cell_type)) {
    stop("Pando multi-task condition fit returned the wrong cell type.",
         call. = FALSE)
  }
  if (!isTRUE(all.equal(as.numeric(fit$padj_threshold), threshold))) {
    stop("Pando returned a condition fit with the wrong BH threshold.",
         call. = FALSE)
  }
  fit$parallel_plan <- list(
    scope = "target",
    target_workers = as.integer(allocated_workers),
    target_parallel = target_parallel,
    target_backend = target_backend,
    bounded_by_regcompass_task_budget = TRUE
  )
  list(cell_type = task$cell_type, grn = fitted, fit = fit)
}

.rc_fit_condition_grns_by_cell_type_budget_impl <-
  .rc_fit_condition_grns_by_cell_type

.rc_fit_condition_grns_by_cell_type <- function(
    object, gem, outdir, pfm = NULL, genome,
    condition_col = "condition", celltype_col = "cell_type",
    cell_type = NULL, rna_assay = "RNA", atac_assay = "ATAC",
    min_cells = 20L,
    pando_initiate_args = list(exclude_exons = TRUE),
    pando_motif_args = list(),
    pando_infer_args = list(
      tf_cor = 0.1, peak_cor = 0.05, adjust_method = "BH",
      padj_threshold = 0.05, rank_action = "mark",
      min_residual_df = 1L
    ),
    save_pando_objects = TRUE, BPPARAM = NULL,
    progress_monitor = NULL,
    species = c("auto", "human", "mouse")) {
  answer <- .rc_fit_condition_grns_by_cell_type_budget_impl(
    object = object, gem = gem, outdir = outdir, pfm = pfm, genome = genome,
    condition_col = condition_col, celltype_col = celltype_col,
    cell_type = cell_type, rna_assay = rna_assay, atac_assay = atac_assay,
    min_cells = min_cells, pando_initiate_args = pando_initiate_args,
    pando_motif_args = pando_motif_args, pando_infer_args = pando_infer_args,
    save_pando_objects = save_pando_objects, BPPARAM = BPPARAM,
    progress_monitor = progress_monitor, species = species
  )
  fits <- answer$condition_grn_fits
  workers_by_type <- if (is.list(fits) && length(fits)) {
    data.frame(
      cell_type = names(fits),
      target_workers = vapply(fits, function(fit) {
        value <- fit$parallel_plan$target_workers %||% 1L
        as.integer(value[[1L]])
      }, integer(1)),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(cell_type = character(), target_workers = integer())
  }
  worker_limit <- if (!is.null(BPPARAM) && !identical(BPPARAM, FALSE)) {
    .rc_bpparam_worker_limit(BPPARAM, default = 1L)
  } else 1L
  old_plan <- answer$pando_execution_summary$parallel_plan %||% list()
  plan <- utils::modifyList(old_plan, list(
    scope = "bounded_hierarchical_cell_type_target",
    workers = worker_limit,
    target_workers_by_cell_type = workers_by_type,
    inner_target_parallel = any(workers_by_type$target_workers > 1L),
    nested_parallel = any(workers_by_type$target_workers > 1L) && nrow(workers_by_type) > 1L,
    nested_backend = if (nrow(workers_by_type) > 1L &&
      any(workers_by_type$target_workers > 1L)) "snow" else "none",
    worker_budget_bounded = TRUE,
    worker_budget_rule =
      "even split of total worker cap across concurrently scheduled condition-GRN cell types"
  ))
  answer$pando_execution_summary$parallel_plan <- plan
  answer$normalization_policy$parallel_contract <- plan
  saveRDS(answer$condition_grn_fits,
          file.path(outdir, "pando_condition_grn_fits.rds"))
  saveRDS(answer, file.path(outdir, "single_cell_grn.rds"))
  answer
}
