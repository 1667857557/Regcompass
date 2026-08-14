.rc_fit_pando_by_celltype_route <- function(
    object, gem, outdir, genome, pfm, species, condition_col, celltype_col,
    condition_types, standard_types, rna_assay, atac_assay,
    extra_args, condition_infer_args, standard_infer_args,
    parallel, BPPARAM, progress_monitor,
    target_rsq_threshold = .RC_PANDO_TARGET_RSQ_THRESHOLD_DEFAULT) {
  if (!length(condition_types) && !length(standard_types)) {
    stop("No Pando cell-type job was selected.", call. = FALSE)
  }
  target_rsq_threshold <- .rc_target_rsq_threshold(target_rsq_threshold)
  base <- list(
    gem = gem, outdir = outdir, genome = genome,
    pfm = pfm, species = species, condition_col = condition_col,
    celltype_col = celltype_col, rna_assay = rna_assay,
    atac_assay = atac_assay
  )
  standard_method <- as.character(standard_infer_args$method %||% "ridge")
  standard_ridge <- identical(standard_method, "ridge")
  worker_limit <- if (isTRUE(parallel) && !is.null(BPPARAM) &&
      !identical(BPPARAM, FALSE)) {
    .rc_bpparam_worker_limit(BPPARAM, default = 1L)
  } else 1L

  .rc_step_monitor_event(
    progress_monitor, "cell_type_execution_plan",
    paste(
      "condition GRNs use pooled/global plus condition exact-union no-fusion ridge;",
      "standard Pando uses", standard_method,
      if (standard_ridge) "through the same ridge solver K=1" else ""
    ),
    current = 5L,
    context = list(
      condition_cell_types = length(condition_types),
      standard_cell_types = length(standard_types),
      condition_parallel_scope = if (length(condition_types) &&
          isTRUE(parallel) && worker_limit > 1L) "target" else "serial",
      standard_parallel_scope = if (length(standard_types) && standard_ridge &&
          isTRUE(parallel) && worker_limit > 1L) {
        "target"
      } else if (length(standard_types) > 1L && !standard_ridge &&
                 isTRUE(parallel)) {
        "cell_type"
      } else {
        "serial"
      },
      target_rsq_threshold = target_rsq_threshold,
      nested_parallel = FALSE,
      memory_policy = paste(
        "single parallel level; one ridge cell type resident at a time;",
        "target-specific Pando worker payloads"
      )
    )
  )

  condition_result <- .rc_run_condition_pando_batch(
    object = object,
    condition_types = condition_types,
    base = base,
    extra_args = extra_args,
    condition_infer_args = condition_infer_args,
    parallel = parallel,
    BPPARAM = BPPARAM,
    progress_monitor = progress_monitor
  )
  invisible(gc(verbose = FALSE, full = TRUE))

  standard_values <- list()
  standard_outer_parallel <- FALSE
  if (length(standard_types)) {
    if (standard_ridge) {
      for (i in seq_along(standard_types)) {
        type <- standard_types[[i]]
        cells <- rownames(object@meta.data)[
          as.character(object@meta.data[[celltype_col]]) == type
        ]
        one <- subset(object, cells = cells)
        one <- .rc_stage1_pando_working_object(
          one, rna_assay = rna_assay, atac_assay = atac_assay
        )
        job <- list(cell_type = type, object = one)
        executed <- .rc_run_standard_pando_celltype_job(
          job = job,
          base = base,
          extra_args = extra_args,
          standard_infer_args = standard_infer_args,
          outer_parallel = FALSE,
          progress_monitor = progress_monitor,
          PANDO_BPPARAM = if (isTRUE(parallel)) BPPARAM else NULL
        )
        standard_values[[type]] <- executed$result
        one <- NULL
        job <- NULL
        executed <- NULL
        invisible(gc(verbose = FALSE, full = TRUE))
      }
    } else {
      standard_source <- .rc_stage1_pando_working_object(
        object, rna_assay = rna_assay, atac_assay = atac_assay
      )
      standard_inputs <- lapply(standard_types, function(type) {
        cells <- rownames(standard_source@meta.data)[
          as.character(standard_source@meta.data[[celltype_col]]) == type
        ]
        list(
          cell_type = type,
          object = subset(standard_source, cells = cells)
        )
      })
      standard_outer_parallel <- isTRUE(parallel) &&
        length(standard_inputs) > 1L
      executed <- rc_parallel_lapply(
        standard_inputs,
        .rc_run_standard_pando_celltype_job,
        BPPARAM = if (standard_outer_parallel) BPPARAM else FALSE,
        base = base,
        extra_args = extra_args,
        standard_infer_args = standard_infer_args,
        outer_parallel = standard_outer_parallel,
        progress_monitor = progress_monitor,
        PANDO_BPPARAM = NULL
      )
      standard_values <- lapply(executed, `[[`, "result")
      names(standard_values) <- vapply(
        executed, `[[`, character(1), "cell_type"
      )
      standard_source <- NULL
      standard_inputs <- NULL
      executed <- NULL
      invisible(gc(verbose = FALSE, full = TRUE))
    }
  }

  answer <- .rc_merge_pando_results(
    condition_result = condition_result,
    standard_results = standard_values,
    condition_types = condition_types,
    standard_types = standard_types,
    condition_col = condition_col,
    celltype_col = celltype_col,
    outdir = outdir,
    target_rsq_threshold = target_rsq_threshold
  )
  condition_plan <- if (!is.null(condition_result) &&
      is.list(condition_result$pando_execution_summary)) {
    condition_result$pando_execution_summary$parallel_plan %||% list()
  } else {
    list()
  }
  answer$pando_execution_plan <- list(
    scope = if (length(condition_types) && length(standard_types)) {
      "condition_then_standard_target_parallel"
    } else if (length(condition_types)) {
      "condition_target_parallel"
    } else {
      "standard_target_parallel"
    },
    condition_cell_types = condition_types,
    standard_cell_types = standard_types,
    condition_parallel_plan = condition_plan,
    standard_method = standard_method,
    standard_outer_parallel = standard_outer_parallel,
    standard_target_parallel = standard_ridge && isTRUE(parallel) &&
      worker_limit > 1L,
    target_rsq_threshold = target_rsq_threshold,
    target_quality_metric = "selected_lambda_final_full_data_rsq",
    oof_rsq_role = "diagnostic_only",
    nested_parallel = FALSE,
    worker_budget_shared_sequentially = TRUE,
    memory_policy = paste(
      "one ridge cell type resident at a time; target-specific Pando payloads;",
      "worker and batch temporaries released after completion"
    )
  )
  answer
}
