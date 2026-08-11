# Bounded condition target routing. A cell-type GRN budget includes the outer
# cell-type worker itself; only the remaining slots may become Pando target
# workers. The canonical condition-GRN orchestrator remains unchanged and calls
# this task helper through the alias at the end of the file.

.rc_condition_multitask_fit_task_bounded <- function(
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

  grn_worker_budget <- suppressWarnings(as.integer(
    Sys.getenv("REGCOMPASS_TASK_WORKERS", unset = "1")
  ))
  if (!is.finite(grn_worker_budget) || grn_worker_budget < 1L) {
    grn_worker_budget <- 1L
  }
  target_param <- NULL
  target_workers <- 1L
  target_backend <- "serial"
  nested_pool <- FALSE

  if (isTRUE(inner_parallel) && !is.null(PANDO_BPPARAM) &&
      !identical(PANDO_BPPARAM, FALSE)) {
    target_param <- PANDO_BPPARAM
    target_workers <- .rc_bpparam_worker_limit(PANDO_BPPARAM, default = 1L)
    grn_worker_budget <- target_workers
    target_backend <- .rc_bpparam_backend(PANDO_BPPARAM)
  } else {
    child_workers <- max(0L, grn_worker_budget - 1L)
    if (child_workers >= 2L) {
      target_param <- .rc_condition_nested_target_bpparam(child_workers)
      target_workers <- child_workers
      target_backend <- "snow"
      nested_pool <- TRUE
    }
  }
  target_parallel <- !is.null(target_param) && target_workers > 1L

  old_r_libs <- Sys.getenv("R_LIBS", unset = NA_character_)
  on.exit({
    if (target_parallel && !is.null(target_param) &&
        requireNamespace("BiocParallel", quietly = TRUE) &&
        isTRUE(tryCatch(BiocParallel::bpisup(target_param),
                        error = function(e) FALSE))) {
      try(BiocParallel::bpstop(target_param), silent = TRUE)
    }
    if (target_parallel) {
      if (is.na(old_r_libs)) Sys.unsetenv("R_LIBS") else Sys.setenv(R_LIBS = old_r_libs)
    }
    invisible(gc(verbose = FALSE, full = TRUE))
  }, add = TRUE)
  if (target_parallel) {
    Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
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
    grn_worker_budget = as.integer(grn_worker_budget),
    target_workers = as.integer(target_workers),
    target_parallel = target_parallel,
    target_backend = target_backend,
    nested_pool = nested_pool,
    outer_worker_included_in_grn_budget = nested_pool,
    target_pool_released_after_cell_type = TRUE,
    bounded_by_regcompass_task_budget = TRUE
  )
  list(cell_type = task$cell_type, grn = fitted, fit = fit)
}

.rc_condition_multitask_fit_task <- .rc_condition_multitask_fit_task_bounded
