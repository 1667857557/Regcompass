# Route Stage 1 Pando arguments and schedule condition/standard Pando work.

.rc_pando_infer_arg_catalog <- function() {
  list(
    shared = c("tf_cor", "peak_cor", "adjust_method"),
    condition = c(
      "padj_threshold", "rank_action", "min_residual_df",
      "rna_layer", "peak_layer", "peak_value_type"
    ),
    standard = c(
      "peak_to_gene_method", "upstream", "downstream", "extend",
      "only_tss", "method", "alpha", "family", "interaction_term",
      "scale", "aggregate_rna_col", "aggregate_peaks_col", "verbose",
      "maxit", "epsilon", "control", "nlambda", "lambda",
      "lambda.min.ratio", "standardize", "nfolds", "type.measure",
      "solver", "bagging_number", "n_jobs", "p_method", "prior",
      "chains", "cores", "iter", "seed", "params", "nrounds",
      "nthread"
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
      "condition_grn_fit", "project_condition_grn_cells",
      "validate_condition_grn_projection_membership",
      "aggregate_condition_grn_projection",
      "discover_grn_edges", "union_grn_edges", "fit_grn_from_edges",
      "GetNetwork", "gof"
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
    tf_cor = 0.1,
    peak_cor = 0.05,
    adjust_method = "BH",
    padj_threshold = 0.05,
    rank_action = "mark",
    min_residual_df = 1L,
    rna_layer = "data",
    peak_layer = "data",
    peak_value_type = "normalized"
  ), condition_args)
  if (length(condition_types) &&
      (!identical(toupper(as.character(condition_args$adjust_method)), "BH") ||
       !isTRUE(all.equal(as.numeric(condition_args$padj_threshold), 0.05)))) {
    stop("Canonical RegCompass condition effects require BH padj < 0.05.",
         call. = FALSE)
  }

  standard_args <- args[intersect(names(args), standard_allowed)]
  standard_args <- utils::modifyList(list(
    tf_cor = 0.1,
    peak_cor = 0.05,
    adjust_method = "BH"
  ), standard_args)

  disabled_condition <- if (length(standard_types)) {
    intersect(names(args), catalog$condition)
  } else character()
  disabled_standard <- if (length(condition_types)) {
    intersect(names(args), catalog$standard)
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

.rc_standard_pando_sample_size_gate <- function(
    args, n_cells, alpha = 0.05) {
  if (!is.list(args)) {
    stop("Standard Pando inference arguments must be a list.", call. = FALSE)
  }
  if (!is.numeric(n_cells) || length(n_cells) != 1L ||
      !is.finite(n_cells) || n_cells < 4) {
    stop("`n_cells` must be one finite value >= 4.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be one finite value in (0, 1).", call. = FALSE)
  }

  requested_tf_cor <- args$tf_cor %||% 0.1
  if (!is.numeric(requested_tf_cor) || length(requested_tf_cor) != 1L ||
      !is.finite(requested_tf_cor) || requested_tf_cor < 0 ||
      requested_tf_cor >= 1) {
    stop("Standard Pando `tf_cor` must be one finite value in [0, 1).",
         call. = FALSE)
  }

  df <- as.numeric(n_cells) - 2
  t_critical <- stats::qt(1 - alpha / 2, df = df)
  sample_size_floor <- t_critical / sqrt(t_critical^2 + df)
  effective_tf_cor <- max(as.numeric(requested_tf_cor), sample_size_floor)
  args$tf_cor <- effective_tf_cor
  attr(args, "sample_size_aware_tf_cor_gate") <- list(
    method = "two_sided_pearson_t_critical",
    alpha = as.numeric(alpha),
    n_cells = as.integer(n_cells),
    requested_tf_cor = as.numeric(requested_tf_cor),
    sample_size_floor = as.numeric(sample_size_floor),
    effective_tf_cor = as.numeric(effective_tf_cor)
  )
  args
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
    outer_parallel, progress_monitor) {
  if (!is.list(job) || !inherits(job$object, "Seurat")) {
    stop("Invalid standard-Pando cell-type job.", call. = FALSE)
  }
  thread_state <- .rc_set_internal_single_thread()
  on.exit(.rc_restore_internal_threads(thread_state), add = TRUE)

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
  gated_infer_args <- .rc_standard_pando_sample_size_gate(
    standard_infer_args, n_cells = ncol(job$object)
  )
  gate <- attr(
    gated_infer_args, "sample_size_aware_tf_cor_gate", exact = TRUE
  )
  attr(gated_infer_args, "sample_size_aware_tf_cor_gate") <- NULL
  gated_infer_args <- .rc_standard_pando_single_process_args(
    gated_infer_args
  )
  pando_parallel_contract <- attr(
    gated_infer_args, "regcompass_pando_parallel_contract", exact = TRUE
  )
  attr(gated_infer_args, "regcompass_pando_parallel_contract") <- NULL
  args$pando_infer_args <- gated_infer_args
  args$parallel <- FALSE
  value <- do.call(.rc_fit_standard_pando_by_cell_type, args)
  if (is.list(value$normalization_policy)) {
    value$normalization_policy$sample_size_aware_tf_cor_gate <- gate
    value$normalization_policy$parallel_contract <-
      pando_parallel_contract
  }
  if (is.data.frame(value$condition_fit_status) &&
      nrow(value$condition_fit_status)) {
    value$condition_fit_status$tf_cor_requested <- gate$requested_tf_cor
    value$condition_fit_status$tf_cor_sample_size_floor <-
      gate$sample_size_floor
    value$condition_fit_status$tf_cor_effective <- gate$effective_tf_cor
    value$condition_fit_status$tf_cor_gate_alpha <- gate$alpha
  }
  list(cell_type = job$cell_type, route = "standard_pando", result = value)
}

.rc_run_condition_pando_batch <- function(
    object, condition_types, base, extra_args, condition_infer_args,
    parallel, BPPARAM, progress_monitor) {
  if (!length(condition_types)) return(NULL)
  cells <- rownames(object@meta.data)[
    as.character(object@meta.data[[base$celltype_col]]) %in% condition_types
  ]
  if (!length(cells)) {
    stop("No cells remain for condition-GRN cell types.", call. = FALSE)
  }
  condition_object <- subset(object, cells = cells)
  args <- c(base[setdiff(names(base), names(extra_args))], extra_args)
  args$object <- condition_object
  args$cell_type <- condition_types
  args$outdir <- file.path(base$outdir, "condition")
  args$pando_infer_args <- condition_infer_args
  args$BPPARAM <- if (isTRUE(parallel)) BPPARAM else FALSE
  args$progress_monitor <- progress_monitor
  do.call(.rc_fit_condition_grns_by_cell_type, args)
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

  .rc_step_monitor_event(
    progress_monitor, "cell_type_execution_plan",
    paste(
      "condition GRNs use condition x cell-type tasks with exact-dictionary",
      "barriers; standard Pando uses broad-cell-type jobs"
    ),
    current = 5L,
    context = list(
      condition_cell_types = length(condition_types),
      standard_cell_types = length(standard_types),
      condition_parallel_scope = if (length(condition_types)) {
        "condition_x_cell_type"
      } else {
        "not_applicable"
      },
      standard_parallel_scope = if (length(standard_types)) {
        "cell_type"
      } else {
        "not_applicable"
      },
      nested_parallel = FALSE
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
      list(
        cell_type = type,
        object = subset(object, cells = cells)
      )
    })
    executed <- rc_parallel_lapply(
      standard_inputs,
      .rc_run_standard_pando_celltype_job,
      BPPARAM = if (standard_outer_parallel) BPPARAM else FALSE,
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
  answer$pando_execution_plan <- list(
    scope = if (length(condition_types) && length(standard_types)) {
      "condition_x_cell_type_then_standard_cell_type"
    } else if (length(condition_types)) {
      "condition_x_cell_type"
    } else {
      "standard_cell_type"
    },
    condition_cell_types = condition_types,
    standard_cell_types = standard_types,
    condition_parallel_plan = condition_plan,
    standard_outer_parallel = standard_outer_parallel,
    nested_parallel = FALSE,
    worker_budget_shared_sequentially = TRUE
  )
  answer
}
