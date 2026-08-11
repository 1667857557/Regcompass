# Optional standard-Pando ridge routing. Original GLM remains the default.

.rc_pando_infer_arg_catalog_standard_ridge_impl <- .rc_pando_infer_arg_catalog
.rc_pando_infer_arg_catalog <- function() {
  out <- .rc_pando_infer_arg_catalog_standard_ridge_impl()
  out$standard <- unique(c(
    out$standard,
    "ridge_control", "rank_action", "min_residual_df"
  ))
  out
}

.rc_route_pando_infer_args_standard_ridge_impl <- .rc_route_pando_infer_args
.rc_route_pando_infer_args <- function(
    args, condition_types = character(), standard_types = character()) {
  if (!is.list(args)) {
    stop("`pando_infer_args` must be a list.", call. = FALSE)
  }
  answer <- .rc_route_pando_infer_args_standard_ridge_impl(
    args,
    condition_types = condition_types,
    standard_types = standard_types
  )
  method <- as.character(answer$standard$method %||% "glm")
  if (length(method) != 1L || is.na(method) || !nzchar(method)) {
    stop("Standard Pando `method` must be one non-empty value.", call. = FALSE)
  }
  ridge_fields <- c("ridge_control", "rank_action", "min_residual_df")
  shared_ridge_fields <- c("rank_action", "min_residual_df")

  if (identical(method, "ridge")) {
    answer$standard$ridge_control <- answer$standard$ridge_control %||% list()
    if (!is.list(answer$standard$ridge_control)) {
      stop("Standard Pando `ridge_control` must be a list.", call. = FALSE)
    }
    answer$standard$rank_action <- answer$standard$rank_action %||% "mark"
    answer$standard$min_residual_df <-
      answer$standard$min_residual_df %||% 1L
  } else {
    if ("ridge_control" %in% names(args) && length(args$ridge_control)) {
      stop(
        "`pando_infer_args$ridge_control` requires standard Pando ",
        "`method = \"ridge\"`.", call. = FALSE
      )
    }
    answer$standard[ridge_fields] <- NULL
  }

  if (is.data.frame(answer$diagnostics) && nrow(answer$diagnostics)) {
    false_condition_disable <-
      answer$diagnostics$route == "condition_grn" &
      answer$diagnostics$argument %in% shared_ridge_fields
    valid_standard_ridge <- identical(method, "ridge") &
      answer$diagnostics$route == "standard_pando" &
      answer$diagnostics$argument %in% ridge_fields
    answer$diagnostics <- answer$diagnostics[
      !(false_condition_disable | valid_standard_ridge),
      , drop = FALSE
    ]
    rownames(answer$diagnostics) <- NULL
  }
  answer
}

.rc_standard_pando_infer_args_ridge_impl <- .rc_standard_pando_infer_args
.rc_standard_pando_infer_args <- function(args) {
  if (!is.list(args)) stop("`pando_infer_args` must be a list.", call. = FALSE)
  method <- as.character(args$method %||% "glm")
  is_ridge <- length(method) == 1L && !is.na(method) && identical(method, "ridge")
  target_param <- args$BPPARAM %||% NULL
  if (is_ridge) args$BPPARAM <- NULL
  answer <- .rc_standard_pando_infer_args_ridge_impl(args)
  if (is_ridge) {
    answer$padj_threshold <- .rc_standard_pando_padj_fixed
  }
  if (is_ridge && !is.null(target_param) && !identical(target_param, FALSE)) {
    if (!requireNamespace("BiocParallel", quietly = TRUE) ||
        !methods::is(target_param, "BiocParallelParam")) {
      stop("Standard ridge `BPPARAM` must be a BiocParallelParam.",
           call. = FALSE)
    }
    answer$BPPARAM <- target_param
  }
  answer
}

.rc_run_standard_pando_celltype_job_glm_impl <-
  .rc_run_standard_pando_celltype_job

.rc_run_standard_pando_celltype_job <- function(
    job, base, extra_args, standard_infer_args,
    outer_parallel, progress_monitor) {
  method <- as.character(standard_infer_args$method %||% "glm")
  if (!identical(method, "ridge")) {
    return(.rc_run_standard_pando_celltype_job_glm_impl(
      job = job,
      base = base,
      extra_args = extra_args,
      standard_infer_args = standard_infer_args,
      outer_parallel = outer_parallel,
      progress_monitor = progress_monitor
    ))
  }
  if (!is.list(job) || !inherits(job$object, "Seurat")) {
    stop("Invalid standard-Pando cell-type job.", call. = FALSE)
  }

  thread_state <- .rc_set_internal_single_thread()
  target_param <- NULL
  on.exit({
    if (!is.null(target_param) &&
        requireNamespace("BiocParallel", quietly = TRUE) &&
        isTRUE(tryCatch(BiocParallel::bpisup(target_param),
                        error = function(e) FALSE))) {
      try(BiocParallel::bpstop(target_param), silent = TRUE)
    }
    .rc_restore_internal_threads(thread_state)
    invisible(gc(verbose = FALSE, full = TRUE))
  }, add = TRUE)

  grn_budget <- suppressWarnings(as.integer(
    Sys.getenv("REGCOMPASS_TASK_WORKERS", unset = "1")
  ))
  if (!is.finite(grn_budget) || grn_budget < 1L) grn_budget <- 1L
  target_workers <- if (isTRUE(outer_parallel)) {
    max(0L, grn_budget - 1L)
  } else {
    grn_budget
  }
  target_parallel <- target_workers >= 2L
  if (target_parallel) {
    target_param <- .rc_condition_nested_target_bpparam(target_workers)
  }

  job_extra <- extra_args
  motif_args <- job_extra$pando_motif_args %||% list()
  if (outer_parallel && is.list(motif_args) &&
      !is.null(motif_args$cache_dir)) {
    motif_args$cache_dir <- file.path(
      motif_args$cache_dir, .rc_safe_path_component(job$cell_type)
    )
    job_extra$pando_motif_args <- motif_args
  }
  args <- c(base[setdiff(names(base), names(job_extra))], job_extra)
  args$object <- job$object
  args$cell_type <- job$cell_type
  args$progress_monitor <- if (outer_parallel) NULL else progress_monitor
  args$outdir <- file.path(
    base$outdir, "standard", .rc_safe_path_component(job$cell_type)
  )

  fixed_infer_args <- .rc_standard_pando_single_process_args(
    standard_infer_args
  )
  original_contract <- attr(
    fixed_infer_args, "regcompass_pando_parallel_contract", exact = TRUE
  )
  attr(fixed_infer_args, "regcompass_pando_parallel_contract") <- NULL
  if (target_parallel) fixed_infer_args$BPPARAM <- target_param
  args$pando_infer_args <- fixed_infer_args
  args$parallel <- target_parallel

  value <- do.call(.rc_fit_standard_pando_by_cell_type, args)
  parallel_contract <- utils::modifyList(
    original_contract %||% list(),
    list(
      scope = "standard_pando_ridge_target",
      method = "ridge",
      grn_worker_budget = as.integer(grn_budget),
      outer_celltype_worker = isTRUE(outer_parallel),
      target_workers = as.integer(if (target_parallel) target_workers else 1L),
      target_parallel = target_parallel,
      target_backend = if (target_parallel) "snow" else "serial",
      target_pool_released_after_cell_type = TRUE,
      worker_budget_bounded = TRUE
    )
  )
  if (is.list(value$normalization_policy)) {
    value$normalization_policy$parallel_contract <- parallel_contract
  }
  list(cell_type = job$cell_type, route = "standard_pando", result = value)
}

.rc_fit_pando_by_celltype_route_standard_ridge_impl <-
  .rc_fit_pando_by_celltype_route

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
  standard_method <- as.character(standard_infer_args$method %||% "glm")
  standard_ridge <- identical(standard_method, "ridge")

  .rc_step_monitor_event(
    progress_monitor, "cell_type_execution_plan",
    paste(
      "condition GRNs use native Pando significant-union multi-task ridge;",
      "standard Pando uses", standard_method, "broad-cell-type jobs"
    ),
    current = 5L,
    context = list(
      condition_cell_types = length(condition_types),
      standard_cell_types = length(standard_types),
      condition_parallel_scope = if (length(condition_types)) {
        "bounded_cell_type_plus_target"
      } else {
        "not_applicable"
      },
      standard_parallel_scope = if (length(standard_types)) {
        if (standard_ridge) "bounded_cell_type_plus_target" else "cell_type"
      } else {
        "not_applicable"
      },
      worker_budget_shared_sequentially = TRUE
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
  standard_outer_parallel <- isTRUE(parallel) && length(standard_types) > 1L
  if (length(standard_types)) {
    standard_inputs <- lapply(standard_types, function(type) {
      cells <- rownames(object@meta.data)[
        as.character(object@meta.data[[celltype_col]]) == type
      ]
      list(cell_type = type, object = subset(object, cells = cells))
    })
    standard_dispatch_param <- if (isTRUE(parallel) && standard_ridge) {
      BPPARAM
    } else if (standard_outer_parallel) {
      BPPARAM
    } else {
      FALSE
    }
    executed <- rc_parallel_lapply(
      standard_inputs,
      .rc_run_standard_pando_celltype_job,
      BPPARAM = standard_dispatch_param,
      base = base,
      extra_args = extra_args,
      standard_infer_args = standard_infer_args,
      outer_parallel = standard_outer_parallel,
      progress_monitor = progress_monitor
    )
    standard_values <- lapply(executed, `[[`, "result")
    names(standard_values) <- vapply(
      executed, `[[`, character(1), "cell_type"
    )
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
  standard_allocation <- if (length(standard_values)) {
    do.call(rbind, lapply(names(standard_values), function(type) {
      contract <- standard_values[[type]]$normalization_policy$parallel_contract %||%
        list()
      data.frame(
        cell_type = type,
        method = standard_method,
        grn_worker_budget = as.integer(contract$grn_worker_budget %||% 1L),
        target_workers = as.integer(contract$target_workers %||% 1L),
        target_parallel = isTRUE(contract$target_parallel),
        stringsAsFactors = FALSE
      )
    }))
  } else {
    data.frame(
      cell_type = character(), method = character(),
      grn_worker_budget = integer(), target_workers = integer(),
      target_parallel = logical()
    )
  }
  answer$pando_execution_plan <- list(
    scope = if (length(condition_types) && length(standard_types)) {
      "condition_multitask_then_standard_cell_type"
    } else if (length(condition_types)) {
      "condition_multitask"
    } else {
      "standard_cell_type"
    },
    condition_cell_types = condition_types,
    standard_cell_types = standard_types,
    condition_parallel_plan = condition_plan,
    standard_method = standard_method,
    standard_outer_parallel = standard_outer_parallel,
    standard_target_parallel = any(standard_allocation$target_parallel),
    standard_worker_allocation = standard_allocation,
    nested_parallel = isTRUE(condition_plan$nested_parallel) ||
      (standard_outer_parallel && any(standard_allocation$target_parallel)),
    worker_budget_shared_sequentially = TRUE,
    target_pool_release_policy = "release_after_each_cell_type_fit"
  )
  answer
}
