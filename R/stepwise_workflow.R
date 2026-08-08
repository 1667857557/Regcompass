.rc_workflow_signature <- function(x) {
  params <- x$params %||% x$workflow_params %||% list()
  params[c(
    "condition_col", "celltype_col", "cell_type", "rna_assay", "atac_assay",
    "analysis_mode"
  )]
}

.rc_validate_grn_metacell_group_coverage <- function(
    grn_result, metacell_meta,
    condition_col = "condition", celltype_col = "cell_type") {
  group_cols <- c(condition_col, celltype_col)
  status <- grn_result$condition_fit_status
  if (!is.data.frame(status) ||
      !all(c(group_cols, "status") %in% colnames(status))) {
    stop("GRN status is incomplete for group coverage validation.",
         call. = FALSE)
  }
  if (!is.data.frame(metacell_meta) ||
      !all(group_cols %in% colnames(metacell_meta))) {
    stop("Metacell metadata are incomplete for group coverage validation.",
         call. = FALSE)
  }
  status$.group_id <- rc_make_stratum_id(status, group_cols)
  metacell_meta$.group_id <- rc_make_stratum_id(metacell_meta, group_cols)
  grn_rows <- split(seq_len(nrow(status)), status$.group_id)
  grn_summary <- do.call(rbind, lapply(grn_rows, function(rows) {
    one <- status[rows, , drop = FALSE]
    values <- one[1L, group_cols, drop = FALSE]
    data.frame(
      values,
      group_id = as.character(one$.group_id[[1L]]),
      grn_status = paste(sort(unique(as.character(one$status))), collapse = ";"),
      n_single_cells = sum(as.numeric(one$n_cells %||% 0), na.rm = TRUE),
      n_active_edges = sum(as.numeric(one$n_active_edges %||% 0), na.rm = TRUE),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))
  metacell_rows <- split(seq_len(nrow(metacell_meta)), metacell_meta$.group_id)
  metacell_summary <- do.call(rbind, lapply(metacell_rows, function(rows) {
    one <- metacell_meta[rows, , drop = FALSE]
    values <- one[1L, group_cols, drop = FALSE]
    data.frame(
      values,
      group_id = as.character(one$.group_id[[1L]]),
      n_metacells = nrow(one),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))
  coverage <- merge(
    grn_summary, metacell_summary,
    by = c(group_cols, "group_id"), all = TRUE, sort = TRUE
  )
  coverage$grn_available <- !is.na(coverage$grn_status) &
    coverage$grn_status == "ok"
  coverage$has_active_pando_evidence <-
    !is.na(coverage$n_active_edges) & coverage$n_active_edges > 0
  coverage$metacells_available <-
    !is.na(coverage$n_metacells) & coverage$n_metacells > 0
  coverage$coverage_complete <- coverage$grn_available &
    coverage$metacells_available
  invalid <- coverage[!coverage$coverage_complete, , drop = FALSE]
  if (nrow(invalid)) {
    stop(
      "GRN and metacell groups do not align: ",
      paste(invalid$group_id, collapse = "; "),
      call. = FALSE
    )
  }
  rownames(coverage) <- NULL
  coverage
}

#' Build condition-pure metacells from the Stage 1 analysis cell set
#'
#' Each broad cell type receives one independent multimodal WNN graph. All
#' conditions within that cell type jointly determine modality weights,
#' neighbours, and Walktrap clusters; condition splits parent membership only
#' after clustering. When a Stage 1 result is supplied, Stage 2 reproduces its
#' exact ordered cell set and validates workflow parameters.
#'
#' When `fragment_files` is supplied, Stage 2 first builds the final SuperCell
#' membership, then aggregates the original single-cell fragments to those
#' metacell barcodes with `SuperCell::AggregateFragmentFile()`. Fragment
#' aggregation uses the same top-level `workers` cap as the rest of RegCompass;
#' `metacell_args$fragment_args$workers` is not a separate public control.
#'
#' @param fragment_files `NULL`/`FALSE` to aggregate the existing ATAC count
#' matrix, a fragment path (or named path vector for multiple samples), or an
#' explicit fragment mapping data frame. Fragment input affects ATAC recounting
#' only; sample identity is not used as a metacell grouping variable.
#' @param workers Total RegCompass worker cap, default 10. The effective cap is
#' `min(workers, max(1, detected logical CPUs - 2))`. Stage 2 uses this value for
#' fragment aggregation and reserves two logical CPUs globally.
#' @param grn Optional output of `rc_regcompass_step_grn()`. Supplying it enables
#' exact Stage 1 cell-ID and parameter validation.
#' @export
rc_regcompass_step_metacells <- function(
    object, outdir,
    condition_col = "condition",
    celltype_col = "cell_type",
    cell_type = NULL,
    rna_assay = "RNA",
    atac_assay = "ATAC",
    fragment_files = NULL,
    metacell_args = list(),
    workers = 10L,
    progress = getOption("RegCompassR.progress", TRUE),
    grn = NULL) {
  monitor <- .rc_step_monitor_start("metacells", outdir, progress)
  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  if (!is.list(metacell_args)) {
    stop("`metacell_args` must be a list.", call. = FALSE)
  }
  worker_config <- .rc_stage_worker_config(workers, argument = "workers")
  worker_limit <- worker_config$worker_limit

  fragment_enabled <- .rc_fragment_input_enabled(fragment_files)
  supplied_fragment_args <- metacell_args$fragment_args %||% list()
  if ("workers" %in% names(supplied_fragment_args)) {
    stop(
      "`metacell_args$fragment_args$workers` has been removed. Use the ",
      "top-level `workers` argument; it is the single RegCompass parallel cap.",
      call. = FALSE
    )
  }
  fragment_args <- .rc_resolve_fragment_aggregation_args(supplied_fragment_args)
  fragment_args$workers <- worker_limit
  metacell_core_args <- metacell_args
  metacell_core_args$fragment_args <- NULL

  n_input <- ncol(object)
  if (is.null(grn)) {
    cell_set <- .rc_build_stage_analysis_cell_set(
      object = object,
      condition_col = condition_col,
      celltype_col = celltype_col,
      cell_type = cell_type,
      pando_args = list(min_cells = .rc_stage1_min_cells_fixed)
    )
    contract <- list(
      source = "independent_stage1_filter_reapplication",
      min_cells = cell_set$min_cells,
      retained_cells = cell_set$retained_cells,
      retained_cell_types = cell_set$retained_cell_types,
      condition_pando_cell_types = cell_set$condition_pando_cell_types,
      standard_pando_cell_types = cell_set$standard_pando_cell_types,
      diagnostics = cell_set$diagnostics,
      analysis_mode = cell_set$analysis_mode,
      condition_levels = cell_set$condition_levels
    )
  } else {
    contract <- .rc_validate_stage1_cell_set(grn)
    expected_condition_col <- if (
      "requested_condition_col" %in% names(grn$params)
    ) {
      grn$params$requested_condition_col
    } else {
      grn$params$condition_col
    }
    expected_celltype_col <- grn$params$celltype_col %||% celltype_col
    expected_rna_assay <- grn$params$rna_assay %||% rna_assay
    expected_atac_assay <- grn$params$atac_assay %||% atac_assay
    if (!identical(condition_col, expected_condition_col)) {
      stop("Stage 2 `condition_col` differs from Stage 1.", call. = FALSE)
    }
    if (!identical(celltype_col, expected_celltype_col)) {
      stop("Stage 2 `celltype_col` differs from Stage 1.", call. = FALSE)
    }
    if (!identical(rna_assay, expected_rna_assay) ||
        !identical(atac_assay, expected_atac_assay)) {
      stop("Stage 2 assay names differ from Stage 1.", call. = FALSE)
    }
    if (!is.null(cell_type) &&
        !setequal(trimws(as.character(cell_type)),
                  contract$retained_cell_types)) {
      stop("Stage 2 `cell_type` differs from the Stage 1 retained types.",
           call. = FALSE)
    }
  }

  object <- .rc_subset_to_stage1_cell_set(object, contract)
  cell_type <- contract$retained_cell_types
  design <- .rc_resolve_condition_design(object, condition_col)
  object <- design$object
  effective_condition_col <- design$condition_col
  resolved_analysis_mode <- if (is.null(grn)) {
    contract$analysis_mode
  } else {
    grn$params$analysis_mode
  }
  resolved_fallback_reason <- if (is.null(grn)) {
    design$fallback_reason
  } else {
    grn$params$fallback_reason
  }
  object <- .rc_prepare_seurat_assays(
    object,
    assays = c(rna_assay, atac_assay),
    required_layers = "counts"
  )
  object <- .rc_clear_signac_fragments(object, atac_assay = atac_assay)

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  pooled <- .rc_make_condition_celltype_metacells(
    object = object,
    outdir = outdir,
    condition_col = effective_condition_col,
    celltype_col = celltype_col,
    cell_type = cell_type,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    metacell_args = metacell_core_args
  )

  fragment_run <- NULL
  if (fragment_enabled) {
    fragment_run <- .rc_aggregate_single_cell_fragments(
      fragment_files = fragment_files,
      membership = pooled$membership,
      outdir = outdir,
      fragment_args = fragment_args
    )
    pooled$metacell_object <- .rc_recount_atac_from_fragment_manifest(
      object = pooled$metacell_object,
      fragment_manifest = fragment_run$fragment_manifest,
      atac_assay = atac_assay,
      require_complete = TRUE,
      process_n = fragment_args$process_n,
      call_peaks = fragment_args$call_peaks,
      macs2_path = fragment_args$macs2_path,
      effective_genome_size = fragment_args$effective_genome_size,
      peak_calling_args = fragment_args$peak_calling_args,
      peak_calling_outdir = file.path(outdir, "fragments", "macs2")
    )
    pooled$fragment_manifest <- fragment_run$fragment_manifest
    pooled$fragment_files <- fragment_run$fragment_files
    pooled$fragment_aggregation <- list(
      schema_version = "regcompass_metacell_fragment_aggregation_v1",
      source = "SuperCell::AggregateFragmentFile",
      args = fragment_run$args,
      details = fragment_run$details
    )
    pooled$input_design$atac_aggregation_method <-
      "single_cell_fragments_to_metacell_then_featurematrix_recount"
    pooled$input_design$fragment_files_supplied <- TRUE
    .rc_write_tsv_gz(
      fragment_run$fragment_manifest,
      file.path(outdir, "fragments", "fragment_manifest.tsv.gz")
    )
  } else {
    pooled$input_design$atac_aggregation_method <-
      "existing_single_cell_ATAC_matrix_sum"
    pooled$input_design$fragment_files_supplied <- FALSE
  }

  metacell_object <- .rc_normalize_condition_metacell_object(
    pooled, rna_assay, atac_assay
  )
  if (!setequal(
    colnames(metacell_object),
    as.character(pooled$metacell_meta$metacell_id)
  )) {
    stop("Metacell object and metadata contain different units.",
         call. = FALSE)
  }

  .rc_write_tsv_gz(
    pooled$metacell_meta, file.path(outdir, "metacell_metadata.tsv.gz")
  )
  .rc_write_tsv_gz(
    pooled$membership, file.path(outdir, "metacell_membership.tsv.gz")
  )
  .rc_write_tsv_gz(
    pooled$celltype_composition,
    file.path(outdir, "metacell_celltype_composition.tsv.gz")
  )
  .rc_write_tsv_gz(
    pooled$celltype_composition_summary,
    file.path(outdir, "metacell_celltype_summary.tsv.gz")
  )

  resolved_metacell_args <- modifyList(
    .rc_condition_metacell_defaults(), metacell_core_args
  )
  resolved_metacell_args$fragment_args <- fragment_args
  answer <- list(
    pooled = pooled,
    metacell_object = metacell_object,
    cell_filter = contract,
    params = list(
      requested_condition_col = design$requested_condition_col,
      condition_col = effective_condition_col,
      condition_levels = design$condition_levels,
      analysis_mode = resolved_analysis_mode,
      fallback_reason = resolved_fallback_reason,
      cell_type_analysis_mode = if (is.null(grn)) NULL else
        grn$params$cell_type_analysis_mode,
      celltype_col = celltype_col,
      cell_type = cell_type,
      rna_assay = rna_assay,
      atac_assay = atac_assay,
      workers = worker_limit,
      parallel_backend = worker_config$actual_backend,
      reserved_cpus = worker_config$reserved_cpus,
      fragment_files_supplied = fragment_enabled,
      fragment_input_type = if (!fragment_enabled) {
        "none"
      } else if (is.data.frame(fragment_files)) {
        "explicit_manifest"
      } else if (length(fragment_files) == 1L) {
        "single_path"
      } else {
        "named_path_vector"
      },
      atac_aggregation_method = pooled$input_design$atac_aggregation_method,
      metacell_args = resolved_metacell_args,
      supercell_api = pooled$input_design$native_supercell_api,
      graph_group_argument = pooled$input_design$graph_group_argument,
      condition_argument = pooled$input_design$condition_argument,
      graph_method = pooled$input_design$graph_method,
      clustering_method = pooled$input_design$clustering_method,
      aggregation_method = pooled$input_design$aggregation_method,
      graph_scope = pooled$input_design$graph_scope,
      condition_scope = pooled$input_design$condition_scope,
      membership_split_timing = pooled$input_design$membership_split_timing,
      modality_weighting = pooled$input_design$modality_weighting,
      temporary_combined_stratum = FALSE,
      seurat_compatibility =
        metacell_object@misc$regcompass_seurat_compatibility,
      n_input_cells = as.integer(n_input),
      n_stage_cells = as.integer(length(contract$retained_cells)),
      cell_set_contract = if (is.null(grn)) {
        "independent_stage1_filter_reapplication_v1"
      } else {
        "stage1_exact_cell_ids_v1"
      }
    )
  )
  answer$cell_filter$stage2_n_input_cells <- as.integer(n_input)
  answer$cell_filter$n_stage_cells <-
    as.integer(length(contract$retained_cells))
  answer$cell_filter$exact_stage1_match <- !is.null(grn)
  class(answer) <- c("regcompass_metacell_step", "list")
  answer <- .rc_step_monitor_finish(answer, monitor)
  saveRDS(answer, file.path(outdir, "step_metacells.rds"))
  answer
}
