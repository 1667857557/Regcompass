# Route Stage 1 Pando arguments and parallelize independent cell-type jobs.

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
    tf_cor = 0.1, peak_cor = 0, adjust_method = "BH",
    padj_threshold = 0.05, rank_action = "mark",
    min_residual_df = 1L
  ), condition_args)
  if (length(condition_types) &&
      (!identical(toupper(as.character(condition_args$adjust_method)), "BH") ||
       !isTRUE(all.equal(as.numeric(condition_args$padj_threshold), 0.05)))) {
    stop("Canonical RegCompass condition effects require BH padj < 0.05.",
         call. = FALSE)
  }

  standard_args <- args[intersect(names(args), standard_allowed)]
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

.rc_merge_pando_grn_data <- function(results, full_object) {
  values <- lapply(results, `[[`, "pando_grn_data")
  values <- values[vapply(values, function(x) inherits(x, "GRNData"), logical(1))]
  if (!length(values)) return(NULL)
  combined <- values[[1L]]
  if (length(values) > 1L) {
    for (value in values[-1L]) {
      for (network_name in names(value@grn@networks)) {
        if (network_name %in% names(combined@grn@networks)) {
          stop("Duplicated Pando network while merging cell-type jobs: ",
               network_name, call. = FALSE)
        }
        combined@grn@networks[[network_name]] <- value@grn@networks[[network_name]]
      }
    }
  }
  if (!inherits(full_object, "Seurat")) {
    stop("The merged condition GRN requires the complete normalized Seurat object.",
         call. = FALSE)
  }
  combined@data <- full_object
  combined@grn@params$condition_grn_fits <- unlist(
    lapply(results, `[[`, "condition_grn_fits"), recursive = FALSE
  )
  combined@grn@params$condition_network_index <- .rc_bind_pando_field(
    results, "pando_network_index"
  )
  network_names <- names(combined@grn@networks)
  if (length(network_names)) combined@grn@active_network <- tail(network_names, 1L)
  combined
}

.rc_merge_condition_job_results <- function(results, full_object) {
  if (!length(results)) return(NULL)
  answer <- results[[1L]]
  if (length(results) > 1L) {
    frame_fields <- c(
      "condition_fit_status", "pando_network_index", "pando_fit_diagnostics",
      "tf_peak_gene_universal", "tf_peak_gene_condition_all",
      "tf_peak_gene_condition", "tf_peak_gene_condition_effect_all",
      "tf_peak_gene_condition_effect", "paired_cell_metadata"
    )
    for (field in frame_fields) {
      answer[[field]] <- .rc_bind_pando_field(results, field)
    }
    if (nrow(answer$paired_cell_metadata)) {
      answer$paired_cell_metadata <- answer$paired_cell_metadata[
        !duplicated(answer$paired_cell_metadata$cell_id), , drop = FALSE
      ]
    }
    answer$paired_cell_ids <- unique(
      as.character(answer$paired_cell_metadata$cell_id)
    )
    answer$target_metabolic_genes <- unique(unlist(
      lapply(results, `[[`, "target_metabolic_genes"), use.names = FALSE
    ))
    answer$condition_grn_fits <- unlist(
      lapply(results, `[[`, "condition_grn_fits"), recursive = FALSE
    )
    summaries <- lapply(results, `[[`, "pando_execution_summary")
    answer$pando_execution_summary <- list(
      fit_engine = paste(unique(unlist(lapply(
        summaries, `[[`, "fit_engine"
      ), use.names = FALSE)), collapse = ";"),
      targets_total = sum(vapply(summaries, function(x) {
        as.integer(x$targets_total %||% 0L)
      }, integer(1))),
      targets_failed = sum(vapply(summaries, function(x) {
        as.integer(x$targets_failed %||% 0L)
      }, integer(1)))
    )
  }
  answer$pando_grn_data <- .rc_merge_pando_grn_data(results, full_object)
  answer
}

.rc_run_pando_celltype_job <- function(
    job, base, extra_args, condition_infer_args, standard_infer_args,
    parallel, outer_parallel, progress_monitor) {
  if (!is.list(job) || !inherits(job$object, "Seurat")) {
    stop("Invalid Pando cell-type job.", call. = FALSE)
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
  if (identical(job$route, "condition_grn")) {
    args$outdir <- file.path(
      base$outdir, "condition", .rc_safe_path_component(job$cell_type)
    )
    args$pando_infer_args <- condition_infer_args
    args$BPPARAM <- FALSE
    value <- do.call(.rc_fit_condition_grns_by_cell_type, args)
  } else if (identical(job$route, "standard_pando")) {
    args$outdir <- file.path(
      base$outdir, "standard", .rc_safe_path_component(job$cell_type)
    )
    args$pando_infer_args <- standard_infer_args
    args$parallel <- isTRUE(parallel) && !outer_parallel
    value <- do.call(.rc_fit_standard_pando_by_cell_type, args)
  } else {
    stop("Unknown Pando route: ", job$route, call. = FALSE)
  }
  list(cell_type = job$cell_type, route = job$route, result = value)
}

.rc_fit_pando_by_celltype_route <- function(
    object, gem, outdir, genome, pfm, species, condition_col, celltype_col,
    condition_types, standard_types, rna_assay, atac_assay,
    extra_args, condition_infer_args, standard_infer_args,
    parallel, BPPARAM, progress_monitor) {
  jobs <- rbind(
    if (length(condition_types)) data.frame(
      cell_type = condition_types, route = "condition_grn",
      stringsAsFactors = FALSE
    ),
    if (length(standard_types)) data.frame(
      cell_type = standard_types, route = "standard_pando",
      stringsAsFactors = FALSE
    )
  )
  if (is.null(jobs) || !nrow(jobs)) {
    stop("No Pando cell-type job was selected.", call. = FALSE)
  }
  outer_parallel <- isTRUE(parallel) && nrow(jobs) > 1L
  .rc_step_monitor_event(
    progress_monitor, "cell_type_execution_plan",
    if (outer_parallel) {
      "parallelizing independent Pando jobs by broad cell type"
    } else {
      "running Pando cell-type jobs without outer parallelism"
    },
    current = 5L,
    context = list(
      jobs = nrow(jobs),
      condition_jobs = sum(jobs$route == "condition_grn"),
      standard_jobs = sum(jobs$route == "standard_pando"),
      outer_parallel = outer_parallel
    )
  )

  job_inputs <- lapply(seq_len(nrow(jobs)), function(index) {
    type <- jobs$cell_type[[index]]
    cells <- rownames(object@meta.data)[
      as.character(object@meta.data[[celltype_col]]) == type
    ]
    list(
      cell_type = type,
      route = jobs$route[[index]],
      object = subset(object, cells = cells)
    )
  })
  base <- list(
    gem = gem, outdir = outdir, genome = genome,
    pfm = pfm, species = species, condition_col = condition_col,
    celltype_col = celltype_col, rna_assay = rna_assay,
    atac_assay = atac_assay
  )
  executed <- rc_parallel_lapply(
    job_inputs,
    .rc_run_pando_celltype_job,
    BPPARAM = if (outer_parallel) BPPARAM else FALSE,
    base = base,
    extra_args = extra_args,
    condition_infer_args = condition_infer_args,
    standard_infer_args = standard_infer_args,
    parallel = parallel,
    outer_parallel = outer_parallel,
    progress_monitor = progress_monitor
  )
  condition_values <- lapply(executed, function(x) {
    if (identical(x$route, "condition_grn")) x$result else NULL
  })
  condition_values <- condition_values[
    !vapply(condition_values, is.null, logical(1))
  ]
  standard_values <- lapply(executed, function(x) {
    if (identical(x$route, "standard_pando")) x$result else NULL
  })
  standard_values <- standard_values[
    !vapply(standard_values, is.null, logical(1))
  ]
  if (length(standard_values)) {
    names(standard_values) <- vapply(executed[
      vapply(executed, function(x) {
        identical(x$route, "standard_pando")
      }, logical(1))
    ], `[[`, character(1), "cell_type")
  }
  condition_result <- .rc_merge_condition_job_results(
    condition_values, full_object = object
  )
  answer <- .rc_merge_pando_results(
    condition_result = condition_result,
    standard_results = standard_values,
    condition_types = condition_types,
    standard_types = standard_types,
    condition_col = condition_col,
    celltype_col = celltype_col,
    outdir = outdir
  )
  single_standard_inner <- isTRUE(parallel) && !outer_parallel &&
    nrow(jobs) == 1L && identical(jobs$route[[1L]], "standard_pando")
  answer$pando_execution_plan <- list(
    scope = if (outer_parallel) {
      "cell_type"
    } else if (single_standard_inner) {
      "target_standard_pando"
    } else {
      "serial"
    },
    n_jobs = nrow(jobs),
    outer_parallel = outer_parallel,
    nested_parallel = FALSE,
    routes = jobs
  )
  answer
}
