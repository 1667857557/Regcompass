.rc_workflow_signature <- function(x) {
  params <- x$params %||% x$workflow_params %||% list()
  params[c("condition_col", "celltype_col", "rna_assay", "atac_assay")]
}

.rc_validate_grn_metacell_group_coverage <- function(
    grn_result, metacell_meta,
    condition_col = "condition", celltype_col = "cell_type") {
  group_cols <- c(condition_col, celltype_col)
  status <- grn_result$group_status
  if (!is.data.frame(status) ||
      !all(c(group_cols, "status") %in% colnames(status))) {
    stop(
      "GRN status is incomplete for condition-by-cell-type coverage validation.",
      call. = FALSE
    )
  }
  if (!is.data.frame(metacell_meta) ||
      !all(group_cols %in% colnames(metacell_meta))) {
    stop(
      "Metacell metadata are incomplete for GRN coverage validation.",
      call. = FALSE
    )
  }
  status$.group_id <- rc_make_stratum_id(status, group_cols)
  metacell_meta$.group_id <- rc_make_stratum_id(metacell_meta, group_cols)

  grn_rows <- split(seq_len(nrow(status)), status$.group_id)
  grn_summary <- do.call(rbind, lapply(grn_rows, function(rows) {
    one <- status[rows, , drop = FALSE]
    values <- one[1L, group_cols, drop = FALSE]
    edge_column <- if ("n_active_edges" %in% colnames(one)) {
      "n_active_edges"
    } else if ("n_significant_edges" %in% colnames(one)) {
      "n_significant_edges"
    } else {
      NULL
    }
    n_edges <- if (is.null(edge_column)) 0 else
      sum(suppressWarnings(as.numeric(one[[edge_column]])), na.rm = TRUE)
    data.frame(
      values,
      group_id = as.character(one$.group_id[[1L]]),
      grn_status = paste(
        sort(unique(as.character(one$status))), collapse = ";"
      ),
      n_single_cells = sum(
        suppressWarnings(as.numeric(one$n_cells %||% 0)), na.rm = TRUE
      ),
      n_active_edges = n_edges,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))

  metacell_rows <- split(
    seq_len(nrow(metacell_meta)), metacell_meta$.group_id
  )
  metacell_summary <- do.call(rbind, lapply(metacell_rows, function(rows) {
    one <- metacell_meta[rows, , drop = FALSE]
    values <- one[1L, group_cols, drop = FALSE]
    purity <- if ("dominant_celltype_fraction" %in% colnames(one)) {
      suppressWarnings(as.numeric(one$dominant_celltype_fraction))
    } else {
      rep(NA_real_, nrow(one))
    }
    mixed <- if ("mixed_celltype_metacell" %in% colnames(one)) {
      one$mixed_celltype_metacell %in% TRUE
    } else {
      rep(FALSE, nrow(one))
    }
    data.frame(
      values,
      group_id = as.character(one$.group_id[[1L]]),
      n_metacells = nrow(one),
      median_dominant_celltype_fraction = if (all(is.na(purity))) {
        NA_real_
      } else {
        stats::median(purity, na.rm = TRUE)
      },
      min_dominant_celltype_fraction = if (all(is.na(purity))) {
        NA_real_
      } else {
        min(purity, na.rm = TRUE)
      },
      n_mixed_celltype_metacells = sum(mixed, na.rm = TRUE),
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
  coverage$has_active_grn_evidence <-
    !is.na(coverage$n_active_edges) & coverage$n_active_edges > 0
  coverage$metacells_available <- !is.na(coverage$n_metacells) &
    coverage$n_metacells > 0
  coverage$coverage_complete <- coverage$grn_available &
    coverage$metacells_available
  invalid <- coverage[!coverage$coverage_complete, , drop = FALSE]
  if (nrow(invalid)) {
    stop(
      "GRN and metacell condition-by-cell-type groups do not align: ",
      paste(invalid$group_id, collapse = "; "),
      paste(
        ". Every scored metacell group requires a successful GRN fit,",
        "and every GRN group requires at least one metacell. A successful",
        "fit may legitimately contain zero active target genes."
      ),
      call. = FALSE
    )
  }
  rownames(coverage) <- NULL
  coverage
}

#' Infer shared-background condition sub-GRNs from single cells
#'
#' The default mode constructs one Pando structural TF-peak-target universe per
#' cell type, then jointly estimates a global edge coefficient and symmetric
#' zero-sum condition deviations with condition-balanced elastic-net regression.
#' Edge stability is estimated by full-size nonparametric bootstrap sampling
#' with replacement within every condition. `legacy_condition_pando` retains
#' the earlier independent condition-by-cell-type Pando fits.
#'
#' @export
rc_regcompass_step_grn <- function(
    object, gem, outdir, genome,
    pfm = NULL,
    species = c("auto", "human", "mouse"),
    condition_col = "condition",
    celltype_col = "cell_type",
    rna_assay = "RNA",
    atac_assay = "ATAC",
    grn_mode = c("multitask_shared_backbone", "legacy_condition_pando"),
    pando_args = list(),
    multitask_args = list(),
    parallel = TRUE,
    BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  monitor <- .rc_step_monitor_start("grn", outdir, progress, total_parts = 5L)
  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)
  grn_mode <- match.arg(grn_mode)
  if (!is.list(pando_args)) stop("`pando_args` must be a list.", call. = FALSE)
  if (!is.list(multitask_args)) {
    stop("`multitask_args` must be a list.", call. = FALSE)
  }
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("`parallel` must be TRUE or FALSE.", call. = FALSE)
  }
  species <- .rc_step_run(monitor, 1L, "validating GEM and species", {
    species <- .rc_infer_gem_species(gem, species)
    rc_validate_gem(gem)
    species
  })
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  object <- .rc_step_run(monitor, 2L, "normalizing RNA and ATAC assays", .rc_normalize_single_cell_grn_object(
    object,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay
  ))
  reserved <- intersect(names(pando_args), c(
    "object", "gem", "outdir", "genome", "pfm", "species",
    "condition_col", "celltype_col", "rna_assay", "atac_assay", "BPPARAM"
  ))
  if (length(reserved)) {
    stop(
      "`pando_args` cannot override workflow fields: ",
      paste(reserved, collapse = ", "), call. = FALSE
    )
  }

  defaults <- list(
    object = object,
    gem = gem,
    outdir = outdir,
    genome = genome,
    pfm = pfm,
    species = species,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )
  if (identical(grn_mode, "multitask_shared_backbone")) {
    legacy_only <- intersect(names(pando_args), c(
      "pando_infer_args", "padj_threshold", "min_abs_estimate",
      "min_model_rsq", "require_padj", "on_group_error"
    ))
    if (length(legacy_only)) {
      stop(
        "The multitask GRN does not accept legacy independent-fit fields: ",
        paste(legacy_only, collapse = ", "),
        ". Use `pando_design_args` and `multitask_args`, or select ",
        "`grn_mode = \"legacy_condition_pando\"`.", call. = FALSE
      )
    }
    defaults$multitask_args <- multitask_args
    defaults$on_celltype_error <- "stop"
    defaults[names(pando_args)] <- NULL
    grn_result <- .rc_step_run(monitor, 3L, "fitting shared multitask GRN", do.call(
      .rc_run_celltype_multitask_grns, c(defaults, pando_args)
    ))
  } else {
    if (length(multitask_args)) {
      warning("`multitask_args` are ignored in legacy GRN mode.", call. = FALSE)
    }
    defaults$on_group_error <- "stop"
    defaults[names(pando_args)] <- NULL
    grn_result <- .rc_step_run(monitor, 3L, "fitting condition-specific GRNs", do.call(
      .rc_run_condition_single_cell_grns, c(defaults, pando_args)
    ))
    if (is.null(grn_result$group_status) &&
        is.data.frame(grn_result$sample_status)) {
      grn_result$group_status <- grn_result$sample_status
      grn_result$sample_status <- NULL
    }
  }
  .rc_step_progress(monitor, 4L, "assembling GRN checkpoint")
  answer <- list(
    grn_result = grn_result,
    gem_fingerprint = .rc_stage_gem_fingerprint(gem),
    params = list(
      condition_col = condition_col,
      celltype_col = celltype_col,
      rna_assay = rna_assay,
      atac_assay = atac_assay,
      grn_mode = grn_mode,
      pando_args = pando_args,
      multitask_args = multitask_args,
      parallel = parallel,
      species = species
    )
  )
  class(answer) <- c("regcompass_grn_step", "list")
  answer <- .rc_step_monitor_finish(answer, monitor)
  .rc_step_progress(monitor, 5L, "saving GRN checkpoint")
  saveRDS(answer, file.path(outdir, "step_grn.rds"))
  answer
}

#' Build condition-only SuperCell2 metacells
#'
#' Cells are stratified only by condition. Cell type is passed as the
#' SuperCell2 label and is retained as metacell provenance; no biological-sample
#' column, balancing rule, or sample-level grouping is used.
#'
#' @export
rc_regcompass_step_metacells <- function(
    object, outdir,
    condition_col = "condition",
    celltype_col = "cell_type",
    rna_assay = "RNA",
    atac_assay = "ATAC",
    fragment_files = FALSE,
    metacell_args = list(),
    progress = getOption("RegCompassR.progress", TRUE)) {
  monitor <- .rc_step_monitor_start("metacells", outdir, progress, total_parts = 5L)
  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)
  object <- .rc_step_run(monitor, 1L, "validating Seurat assay layers", .rc_prepare_seurat_assays(
    object,
    assays = c(rna_assay, atac_assay),
    required_layers = "counts"
  ))
  if (identical(fragment_files, FALSE) || is.null(fragment_files)) {
    object <- .rc_clear_signac_fragments(object, atac_assay = atac_assay)
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  pooled <- .rc_step_run(monitor, 2L, "constructing condition-pooled metacells", .rc_make_condition_pooled_metacells(
    object = object,
    outdir = outdir,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    fragment_files = fragment_files,
    metacell_args = metacell_args
  ))
  metacell_object <- .rc_step_run(monitor, 3L, "normalizing merged metacell assays", .rc_normalize_condition_metacell_object(
    pooled, rna_assay, atac_assay
  ))
  if (!setequal(
    colnames(metacell_object),
    as.character(pooled$metacell_meta$metacell_id)
  )) {
    stop(
      "Merged metacell object and metadata contain different units.",
      call. = FALSE
    )
  }
  .rc_step_progress(monitor, 4L, "writing metacell tables")
  .rc_write_tsv_gz(
    pooled$metacell_meta,
    file.path(outdir, "metacell_metadata.tsv.gz")
  )
  .rc_write_tsv_gz(
    pooled$membership,
    file.path(outdir, "metacell_membership.tsv.gz")
  )
  .rc_write_tsv_gz(
    pooled$celltype_composition,
    file.path(outdir, "metacell_celltype_composition.tsv.gz")
  )
  .rc_write_tsv_gz(
    pooled$celltype_composition_summary,
    file.path(outdir, "metacell_celltype_summary.tsv.gz")
  )
  saveRDS(
    metacell_object,
    file.path(outdir, "merged_metacell_object.rds")
  )
  answer <- list(
    pooled = pooled,
    metacell_object = metacell_object,
    params = list(
      condition_col = condition_col,
      celltype_col = celltype_col,
      rna_assay = rna_assay,
      atac_assay = atac_assay,
      fragment_files = fragment_files,
      metacell_args = modifyList(list(gamma = 30L), metacell_args),
      seurat_compatibility =
        metacell_object@misc$regcompass_seurat_compatibility
    )
  )
  class(answer) <- c("regcompass_metacell_step", "list")
  answer <- .rc_step_monitor_finish(answer, monitor)
  .rc_step_progress(monitor, 5L, "saving metacell checkpoint")
  saveRDS(answer, file.path(outdir, "step_metacells.rds"))
  answer
}
