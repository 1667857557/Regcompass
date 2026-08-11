# Route Stage 1 Pando arguments and schedule condition/standard Pando work.

.rc_stage1_pando_working_object <- function(
    object, rna_assay = "RNA", atac_assay = "ATAC") {
  if (!inherits(object, "Seurat")) {
    stop("Stage 1 Pando working data must inherit from Seurat.", call. = FALSE)
  }
  assays <- unique(c(as.character(rna_assay), as.character(atac_assay)))
  missing <- setdiff(assays, .rc_seurat_assay_names(object))
  if (length(missing)) {
    stop(
      "Stage 1 Pando working data are missing assays: ",
      paste(missing, collapse = ", "), ".", call. = FALSE
    )
  }
  SeuratObject::DefaultAssay(object) <- rna_assay
  slim <- Seurat::DietSeurat(
    object = object,
    counts = TRUE,
    data = TRUE,
    scale.data = FALSE,
    assays = assays,
    dimreducs = NULL,
    graphs = NULL,
    misc = FALSE
  )
  SeuratObject::DefaultAssay(slim) <- rna_assay
  slim
}

.rc_pando_infer_arg_catalog <- function() {
  list(
    shared = c("tf_cor", "peak_cor", "adjust_method"),
    condition = c(
      "padj_threshold", "rank_action", "min_residual_df",
      "rna_layer", "peak_layer", "peak_value_type",
      "condition_ridge_control"
    ),
    standard = c(
      "peak_to_gene_method", "upstream", "downstream", "extend",
      "only_tss", "method", "alpha", "family", "interaction_term",
      "scale", "aggregate_rna_col", "aggregate_peaks_col", "verbose",
      "maxit", "epsilon", "control", "nlambda", "lambda",
      "lambda.min.ratio", "standardize", "nfolds", "type.measure",
      "solver", "bagging_number", "n_jobs", "p_method", "prior",
      "chains", "cores", "iter", "seed", "params", "nrounds",
      "nthread", "ridge_control", "rank_action", "min_residual_df"
    )
  )
}

.rc_route_pando_infer_args <- function(
    args, condition_types = character(), standard_types = character()) {
  if (!is.list(args)) {
    stop("`pando_infer_args` must be a list.", call. = FALSE)
  }
  if (is.null(names(args))) names(args) <- rep("", length(args))
  if (any(!nzchar(names(args)))) {
    stop("Every `pando_infer_args` entry must be named.", call. = FALSE)
  }

  if (length(condition_types)) {
    if (!requireNamespace("Pando", quietly = TRUE)) {
      stop("Package 'Pando' is required for condition-GRN Stage 1.",
           call. = FALSE)
    }
    required_api <- c(
      "infer_condition_grn", "condition_grn_fit", "initiate_grn",
      "find_motifs", "LayerData"
    )
    missing_api <- setdiff(required_api, getNamespaceExports("Pando"))
    if (length(missing_api)) {
      stop(
        "Installed Pando lacks required RegCompass condition-GRN API(s): ",
        paste(missing_api, collapse = ", "), ".", call. = FALSE
      )
    }
  }

  canonical_layers <- list(
    rna_layer = "data",
    peak_layer = "data",
    peak_value_type = "normalized"
  )
  supplied_layers <- intersect(names(args), names(canonical_layers))
  if (length(condition_types) && length(supplied_layers)) {
    invalid <- vapply(supplied_layers, function(field) {
      value <- args[[field]]
      !is.character(value) || length(value) != 1L || is.na(value) ||
        !identical(as.character(value), canonical_layers[[field]])
    }, logical(1))
    if (any(invalid)) {
      stop(
        "Condition-GRN Stage 1 requires rna_layer='data', peak_layer='data', ",
        "and peak_value_type='normalized'. Unsupported override(s): ",
        paste(supplied_layers[invalid], collapse = ", "), ".",
        call. = FALSE
      )
    }
  }

  catalog <- .rc_pando_infer_arg_catalog()
  condition_allowed <- unique(c(catalog$shared, catalog$condition))
  standard_allowed <- unique(c(catalog$shared, catalog$standard))
  known <- unique(c(condition_allowed, standard_allowed))
  unknown <- setdiff(names(args), known)
  if (length(unknown)) {
    stop(
      "Unsupported `pando_infer_args`: ", paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }

  condition_args <- args[intersect(names(args), condition_allowed)]
  condition_args <- utils::modifyList(list(
    tf_cor = 0.05,
    peak_cor = 0.05,
    adjust_method = "BH",
    padj_threshold = 0.05,
    rank_action = "mark",
    min_residual_df = 1L,
    rna_layer = "data",
    peak_layer = "data",
    peak_value_type = "normalized",
    condition_ridge_control = list()
  ), condition_args)
  condition_threshold <- suppressWarnings(as.numeric(
    condition_args$padj_threshold
  ))
  if (length(condition_types) &&
      (!identical(toupper(as.character(condition_args$adjust_method)), "BH") ||
       length(condition_threshold) != 1L || !is.finite(condition_threshold) ||
       condition_threshold <= 0 || condition_threshold >= 1 ||
       !is.list(condition_args$condition_ridge_control))) {
    stop(
      "Canonical RegCompass condition fits require BH adjustment, ",
      "padj_threshold in (0, 1), and condition_ridge_control as a list.",
      call. = FALSE
    )
  }
  condition_args$padj_threshold <- condition_threshold

  standard_args <- args[intersect(names(args), standard_allowed)]
  standard_args <- utils::modifyList(list(
    tf_cor = 0.05,
    peak_cor = 0.05,
    adjust_method = "BH"
  ), standard_args)
  standard_method <- as.character(standard_args$method %||% "ridge")
  if (length(standard_method) != 1L || is.na(standard_method) ||
      !nzchar(standard_method)) {
    stop("Standard Pando `method` must be one non-empty value.", call. = FALSE)
  }
  standard_args$method <- standard_method
  if (identical(standard_method, "ridge")) {
    standard_args$ridge_control <- standard_args$ridge_control %||% list()
    if (!is.list(standard_args$ridge_control)) {
      stop("Standard Pando `ridge_control` must be a list.", call. = FALSE)
    }
    standard_args$rank_action <- standard_args$rank_action %||% "mark"
    standard_args$min_residual_df <- standard_args$min_residual_df %||% 1L
  } else {
    if ("ridge_control" %in% names(args) && length(args$ridge_control)) {
      stop(
        "`pando_infer_args$ridge_control` requires standard Pando ",
        "`method = \"ridge\"`.", call. = FALSE
      )
    }
    standard_args[c("ridge_control", "rank_action", "min_residual_df")] <- NULL
  }

  condition_only <- setdiff(catalog$condition, catalog$standard)
  standard_only <- setdiff(catalog$standard, catalog$condition)
  disabled_condition <- if (length(standard_types)) {
    intersect(names(args), condition_only)
  } else character()
  disabled_standard <- if (length(condition_types)) {
    intersect(names(args), standard_only)
  } else character()
  diagnostics <- rbind(
    if (length(disabled_condition)) data.frame(
      route = "standard_pando", argument = disabled_condition,
      action = "disabled_condition_only", stringsAsFactors = FALSE
    ),
    if (length(disabled_standard)) data.frame(
      route = "condition_grn", argument = disabled_standard,
      action = "disabled_standard_only", stringsAsFactors = FALSE
    )
  )
  if (is.null(diagnostics)) diagnostics <- data.frame(
    route = character(), argument = character(), action = character(),
    stringsAsFactors = FALSE
  )
  if (nrow(diagnostics)) {
    message(
      "Stage 1 routed `pando_infer_args` by analysis mode; disabled: ",
      paste0(diagnostics$route, "{", diagnostics$argument, "}",
             collapse = ", ")
    )
  }
  list(
    condition = condition_args,
    standard = standard_args,
    diagnostics = diagnostics
  )
}

.rc_standard_pando_single_process_args <- function(args) {
  if (!is.list(args)) return(args)
  controls <- intersect(names(args), c("n_jobs", "cores", "nthread"))
  requested <- if (length(controls)) args[controls] else list()
  for (name in controls) args[[name]] <- 1L
  if (is.list(args$params) && "nthread" %in% names(args$params)) {
    requested$params_nthread <- args$params$nthread
    args$params$nthread <- 1L
  }
  attr(args, "regcompass_pando_parallel_contract") <- list(
    scope = "standard_pando_broad_cell_type_jobs",
    infer_grn_parallel = FALSE,
    inner_worker_limit = 1L,
    overridden_controls = requested
  )
  args
}

.rc_run_standard_pando_celltype_job <- function(
    job, base, extra_args, standard_infer_args,
    outer_parallel, progress_monitor, PANDO_BPPARAM = NULL) {
  if (!is.list(job) || !inherits(job$object, "Seurat")) {
    stop("Invalid standard-Pando cell-type job.", call. = FALSE)
  }
  thread_state <- .rc_set_internal_single_thread()
  on.exit({
    .rc_restore_internal_threads(thread_state)
    invisible(gc(verbose = FALSE, full = TRUE))
  }, add = TRUE)

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
  pando_parallel_contract <- attr(
    fixed_infer_args, "regcompass_pando_parallel_contract", exact = TRUE
  )
  attr(fixed_infer_args, "regcompass_pando_parallel_contract") <- NULL
  method <- as.character(fixed_infer_args$method %||% "ridge")
  ridge_inner_parallel <- identical(method, "ridge") && !isTRUE(outer_parallel) &&
    !is.null(PANDO_BPPARAM) && !identical(PANDO_BPPARAM, FALSE) &&
    .rc_bpparam_worker_limit(PANDO_BPPARAM, default = 1L) > 1L
  args$pando_infer_args <- fixed_infer_args
  args$parallel <- ridge_inner_parallel
  args$BPPARAM <- if (ridge_inner_parallel) PANDO_BPPARAM else NULL
  if (is.list(pando_parallel_contract)) {
    pando_parallel_contract$scope <- if (ridge_inner_parallel) {
      "standard_pando_single_celltype_target_parallel"
    } else {
      "standard_pando_celltype_parallel_or_serial"
    }
    pando_parallel_contract$infer_grn_parallel <- ridge_inner_parallel
    pando_parallel_contract$inner_worker_limit <- if (ridge_inner_parallel) {
      .rc_bpparam_worker_limit(PANDO_BPPARAM, default = 1L)
    } else 1L
    pando_parallel_contract$nested_parallel <- FALSE
  }
  value <- do.call(.rc_fit_standard_pando_by_cell_type, args)
  if (is.list(value$normalization_policy)) {
    value$normalization_policy$parallel_contract <- pando_parallel_contract
  }
  list(cell_type = job$cell_type, route = "standard_pando", result = value)
}

.rc_run_condition_pando_batch <- function(
    object, condition_types, base, extra_args, condition_infer_args,
    parallel, BPPARAM, progress_monitor) {
  if (!length(condition_types)) return(NULL)
  worker_limit <- if (isTRUE(parallel) && !is.null(BPPARAM) &&
      !identical(BPPARAM, FALSE)) {
    .rc_bpparam_worker_limit(BPPARAM, default = 1L)
  } else 1L
  values <- vector("list", length(condition_types))
  names(values) <- condition_types
  for (i in seq_along(condition_types)) {
    type <- condition_types[[i]]
    cells <- rownames(object@meta.data)[
      as.character(object@meta.data[[base$celltype_col]]) == type
    ]
    if (!length(cells)) {
      stop("No cells remain for condition-GRN cell type `", type, "`.",
           call. = FALSE)
    }
    one <- subset(object, cells = cells)
    one <- .rc_stage1_pando_working_object(
      one, rna_assay = base$rna_assay, atac_assay = base$atac_assay
    )
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
