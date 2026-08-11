# Workflow policy: standard Pando defaults to the condition-ridge K=1 solver and
# Pando target-level parallelism is preferred over cell-type-level nesting.
#
# Each cell type is kept resident only while it is being fitted. The paired
# Pando runtime supplies target-specific RNA/ATAC/motif/dictionary payloads to
# workers and releases worker/batch temporaries after completion.

.rc_route_pando_infer_args_before_default_ridge <- .rc_route_pando_infer_args
.rc_standard_pando_infer_args_before_default_ridge <- .rc_standard_pando_infer_args

.rc_route_pando_infer_args <- function(
    args, condition_types = character(), standard_types = character()) {
  if (!is.list(args)) {
    stop("`pando_infer_args` must be a list.", call. = FALSE)
  }
  if (length(standard_types) && is.null(args$method)) {
    args$method <- "ridge"
  }
  .rc_route_pando_infer_args_before_default_ridge(
    args = args,
    condition_types = condition_types,
    standard_types = standard_types
  )
}

.rc_standard_pando_infer_args <- function(args) {
  if (!is.list(args)) {
    stop("`pando_infer_args` must be a list.", call. = FALSE)
  }
  if (is.null(args$method)) args$method <- "ridge"
  .rc_standard_pando_infer_args_before_default_ridge(args)
}

.rc_run_condition_pando_batch <- function(
    object, condition_types, base, extra_args, condition_infer_args,
    parallel, BPPARAM, progress_monitor) {
  if (!length(condition_types)) return(NULL)
  working_object <- .rc_stage1_pando_working_object(
    object, rna_assay = base$rna_assay, atac_assay = base$atac_assay
  )
  worker_limit <- if (isTRUE(parallel) && !is.null(BPPARAM) &&
      !identical(BPPARAM, FALSE)) {
    .rc_bpparam_worker_limit(BPPARAM, default = 1L)
  } else 1L
  values <- vector("list", length(condition_types))
  names(values) <- condition_types
  for (i in seq_along(condition_types)) {
    type <- condition_types[[i]]
    cells <- rownames(working_object@meta.data)[
      as.character(working_object@meta.data[[base$celltype_col]]) == type
    ]
    if (!length(cells)) {
      stop("No cells remain for condition-GRN cell type `", type, "`.",
           call. = FALSE)
    }
    one <- subset(working_object, cells = cells)
    args <- c(base[setdiff(names(base), names(extra_args))], extra_args)
    args$object <- one
    args$cell_type <- type
    args$outdir <- file.path(
      base$outdir, "condition", .rc_safe_path_component(type)
    )
    args$pando_infer_args <- condition_infer_args
    args$BPPARAM <- if (isTRUE(parallel)) BPPARAM else FALSE
    args$progress_monitor <- progress_monitor
    values[[i]] <- do.call(.rc_fit_condition_grns_by_cell_type, args)
    one <- NULL
    args <- NULL
    invisible(gc(verbose = FALSE, full = TRUE))
  }
  working_object <- NULL
  invisible(gc(verbose = FALSE, full = TRUE))
  answer <- if (length(values) == 1L) {
    values[[1L]]
  } else {
    .rc_merge_condition_job_results(values)
  }
  plan <- list(
    scope = "cell_type_serial_target_parallel",
    cell_types = condition_types,
    workers = worker_limit,
    outer_celltype_parallel = FALSE,
    inner_target_parallel = isTRUE(parallel) && worker_limit > 1L,
    nested_parallel = FALSE,
    worker_budget_shared_sequentially = TRUE,
    memory_policy = paste(
      "one cell type resident at a time; Pando target workers receive",
      "target-specific compact payloads"
    )
  )
  if (!is.list(answer$pando_execution_summary)) {
    answer$pando_execution_summary <- list()
  }
  answer$pando_execution_summary$parallel_plan <- plan
  if (!is.list(answer$normalization_policy)) {
    answer$normalization_policy <- list()
  }
  answer$normalization_policy$parallel_contract <- plan
  answer
}

.rc_fit_pando_by_celltype_route <- function(
    object, gem, outdir, genome, pfm, species, condition_col, celltype_col,
    condition_types, standard_types, rna_assay, atac_assay,
    extra_args, condition_infer_args, standard_infer_args,
    parallel, BPPARAM, progress_monitor) {
  if (!length(condition_types) && !length(standard_types)) {
    stop("No Pando cell-type job was selected.", call. = FALSE)
  }
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
      "condition GRNs use native Pando significant-union multi-task ridge;",
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
    standard_source <- .rc_stage1_pando_working_object(
      object, rna_assay = rna_assay, atac_assay = atac_assay
    )
    if (standard_ridge) {
      for (i in seq_along(standard_types)) {
        type <- standard_types[[i]]
        cells <- rownames(standard_source@meta.data)[
          as.character(standard_source@meta.data[[celltype_col]]) == type
        ]
        job <- list(
          cell_type = type,
          object = subset(standard_source, cells = cells)
        )
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
        job <- NULL
        executed <- NULL
        invisible(gc(verbose = FALSE, full = TRUE))
      }
    } else {
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
      standard_inputs <- NULL
      executed <- NULL
      invisible(gc(verbose = FALSE, full = TRUE))
    }
    standard_source <- NULL
    invisible(gc(verbose = FALSE, full = TRUE))
  }

  answer <- .rc_merge_pando_results(
    condition_result = condition_result,
    standard_results = standard_values,
    condition_types = condition_types,
    standard_types = standard_types,
    condition_col = condition_col,
    celltype_col = celltype_col,
    outdir = outdir
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
    nested_parallel = FALSE,
    worker_budget_shared_sequentially = TRUE,
    memory_policy = paste(
      "one ridge cell type resident at a time; target-specific Pando payloads;",
      "worker and batch temporaries released after completion"
    )
  )
  answer
}
