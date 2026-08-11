# Optional standard-Pando ridge routing. Original GLM remains the default.
# Canonical direct definitions remain in the original files. This adapter is
# self-contained and installs the extended helpers only through final aliases;
# no previous implementation is saved or wrapped.

.rc_pando_infer_arg_catalog_standard_ridge <- function() {
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
      "scale", "verbose", "maxit", "epsilon", "control", "nlambda",
      "lambda", "lambda.min.ratio", "standardize", "nfolds",
      "type.measure", "solver", "bagging_number", "n_jobs", "p_method",
      "prior", "chains", "cores", "iter", "seed", "params", "nrounds",
      "nthread", "ridge_control", "rank_action", "min_residual_df"
    )
  )
}

.rc_route_pando_infer_args_standard_ridge <- function(
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

  catalog <- .rc_pando_infer_arg_catalog_standard_ridge()
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
    tf_cor = 0.1,
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
    tf_cor = 0.1,
    peak_cor = 0.05,
    adjust_method = "BH"
  ), standard_args)
  standard_method <- as.character(standard_args$method %||% "glm")
  if (length(standard_method) != 1L || is.na(standard_method) ||
      !nzchar(standard_method)) {
    stop("Standard Pando `method` must be one non-empty value.", call. = FALSE)
  }
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

.rc_standard_pando_infer_args_standard_ridge <- function(args) {
  if (!is.list(args)) stop("`pando_infer_args` must be a list.", call. = FALSE)
  method <- as.character(args$method %||% "glm")
  is_ridge <- length(method) == 1L && !is.na(method) && identical(method, "ridge")
  forbidden <- intersect(names(args), c(
    "object", "genes", "network_name", "aggregate_rna_col",
    "aggregate_peaks_col", "parallel"
  ))
  if (length(forbidden)) {
    stop("Standard Pando inference arguments cannot override managed fields: ",
         paste(forbidden, collapse = ", "), ".", call. = FALSE)
  }
  condition_only <- intersect(names(args), c(
    "candidate_screen", "condition_mix", "condition_weight",
    "reference_condition", "comparison_conditions",
    "lambda_min_ratio", "outer_nfolds", "inner_nfolds",
    "lambda_selection", "active_tol", "max_iter",
    "tol_objective", "tol_coef", "BPPARAM"
  ))
  requested_scale <- args$scale %||% NULL
  if (length(condition_only)) args[condition_only] <- NULL
  if (!is.null(args$interaction_term) &&
      !identical(as.character(args$interaction_term), ":")) {
    stop("Standard RegCompass projection requires `interaction_term = ':'`.",
         call. = FALSE)
  }
  if (!is_ridge) {
    args[c("ridge_control", "rank_action", "min_residual_df",
           "padj_threshold")] <- NULL
  } else {
    args$ridge_control <- args$ridge_control %||% list()
    args$rank_action <- args$rank_action %||% "mark"
    args$min_residual_df <- args$min_residual_df %||% 1L
    args$padj_threshold <- .rc_standard_pando_padj_fixed
  }
  args$scale <- FALSE
  answer <- modifyList(list(interaction_term = ":", scale = FALSE), args)
  attr(answer, "standard_fallback_adjustments") <- list(
    dropped_condition_arguments = condition_only,
    scale_forced_false = !identical(requested_scale, FALSE),
    reason = "one_effective_condition_uses_original_pando_projection_scale"
  )
  answer
}

.rc_run_standard_pando_celltype_job_ridge <- function(
    job, base, extra_args, standard_infer_args,
    outer_parallel, progress_monitor) {
  if (!is.list(job) || !inherits(job$object, "Seurat")) {
    stop("Invalid standard-Pando cell-type job.", call. = FALSE)
  }
  method <- as.character(standard_infer_args$method %||% "glm")
  is_ridge <- identical(method, "ridge")
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

  if (is_ridge) {
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
      fixed_infer_args$BPPARAM <- target_param
    }
    args$parallel <- target_parallel
    pando_parallel_contract <- utils::modifyList(
      pando_parallel_contract %||% list(),
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
  } else {
    args$parallel <- FALSE
  }

  args$pando_infer_args <- fixed_infer_args
  value <- do.call(.rc_fit_standard_pando_by_cell_type, args)
  if (is.list(value$normalization_policy)) {
    value$normalization_policy$parallel_contract <- pando_parallel_contract
  }
  list(cell_type = job$cell_type, route = "standard_pando", result = value)
}

.rc_fit_pando_by_celltype_route_ridge <- function(
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
  condition_fits <- if (!is.null(condition_result)) {
    condition_result$condition_grn_fits %||% list()
  } else {
    list()
  }
  if (is.list(condition_fits) && length(condition_fits)) {
    condition_allocation <- data.frame(
      cell_type = names(condition_fits),
      grn_worker_budget = vapply(condition_fits, function(fit) {
        as.integer((fit$parallel_plan$grn_worker_budget %||% 1L)[[1L]])
      }, integer(1)),
      target_workers = vapply(condition_fits, function(fit) {
        as.integer((fit$parallel_plan$target_workers %||% 1L)[[1L]])
      }, integer(1)),
      nested_pool = vapply(condition_fits, function(fit) {
        isTRUE(fit$parallel_plan$nested_pool)
      }, logical(1)),
      stringsAsFactors = FALSE
    )
    condition_plan$scope <- "bounded_hierarchical_cell_type_target"
    condition_plan$worker_allocation_by_cell_type <- condition_allocation
    condition_plan$inner_target_parallel <-
      any(condition_allocation$target_workers > 1L)
    condition_plan$nested_parallel <- any(condition_allocation$nested_pool)
    condition_plan$nested_backend <-
      if (any(condition_allocation$nested_pool)) "snow" else "none"
    condition_plan$worker_budget_bounded <- TRUE
    condition_plan$target_pools_release_policy <-
      "release_after_each_cell_type_fit"
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

.rc_pando_infer_arg_catalog <- .rc_pando_infer_arg_catalog_standard_ridge
.rc_route_pando_infer_args <- .rc_route_pando_infer_args_standard_ridge
.rc_standard_pando_infer_args <- .rc_standard_pando_infer_args_standard_ridge
.rc_run_standard_pando_celltype_job <- .rc_run_standard_pando_celltype_job_ridge
.rc_fit_pando_by_celltype_route <- .rc_fit_pando_by_celltype_route_ridge
