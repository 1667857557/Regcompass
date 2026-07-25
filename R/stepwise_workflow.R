.rc_workflow_signature <- function(x) {
  params <- x$params %||% x$workflow_params %||% list()
  params[c("condition_col", "celltype_col", "rna_assay", "atac_assay")]
}

.rc_validate_grn_metacell_group_coverage <- function(
    grn_result, metacell_meta,
    condition_col = "condition", celltype_col = "cell_type") {
  group_cols <- c(condition_col, celltype_col)
  status <- grn_result$sample_status
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
    data.frame(
      values,
      group_id = as.character(one$.group_id[[1L]]),
      grn_status = paste(
        sort(unique(as.character(one$status))), collapse = ";"
      ),
      n_single_cells = sum(as.numeric(one$n_cells %||% 0), na.rm = TRUE),
      n_significant_edges = sum(
        as.numeric(one$n_significant_edges %||% 0), na.rm = TRUE
      ),
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
    coverage$grn_status == "ok" & coverage$n_significant_edges > 0
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
        ". Every scored metacell group requires a successful GRN with",
        "significant edges, and every GRN group requires at least one metacell."
      ),
      call. = FALSE
    )
  }
  rownames(coverage) <- NULL
  coverage
}

#' Infer condition-by-cell-type Pando GRNs from single cells
#' @export
rc_regcompass_step_grn <- function(
    object, gem, outdir, pfm, genome,
    condition_col = "condition",
    celltype_col = "cell_type",
    rna_assay = "RNA",
    atac_assay = "ATAC",
    pando_args = list(),
    parallel = TRUE,
    BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  monitor <- .rc_step_monitor_start("grn", outdir, progress)
  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)
  if (!is.list(pando_args)) {
    stop("`pando_args` must be a list.", call. = FALSE)
  }
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("`parallel` must be TRUE or FALSE.", call. = FALSE)
  }
  rc_validate_gem(gem)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  object <- .rc_normalize_single_cell_grn_object(
    object,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay
  )
  reserved <- intersect(names(pando_args), c(
    "object", "gem", "outdir", "pfm", "genome", "condition_col",
    "celltype_col", "rna_assay", "atac_assay", "BPPARAM"
  ))
  if (length(reserved)) {
    stop(
      "`pando_args` cannot override workflow fields: ",
      paste(reserved, collapse = ", "),
      call. = FALSE
    )
  }
  defaults <- list(
    object = object,
    gem = gem,
    outdir = outdir,
    pfm = pfm,
    genome = genome,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE,
    on_group_error = "stop"
  )
  defaults[names(pando_args)] <- NULL
  grn_result <- do.call(
    .rc_run_condition_single_cell_grns, c(defaults, pando_args)
  )
  answer <- list(
    grn_result = grn_result,
    gem_fingerprint = .rc_stage_gem_fingerprint(gem),
    params = list(
      condition_col = condition_col,
      celltype_col = celltype_col,
      rna_assay = rna_assay,
      atac_assay = atac_assay,
      pando_args = pando_args,
      parallel = parallel
    )
  )
  class(answer) <- c("regcompass_grn_step", "list")
  answer <- .rc_step_monitor_finish(answer, monitor)
  saveRDS(answer, file.path(outdir, "step_grn.rds"))
  answer
}

#' Build condition-only SuperCell2 metacells
#' @export
rc_regcompass_step_metacells <- function(
    object, outdir,
    sample_col = NULL,
    condition_col = "condition",
    celltype_col = "cell_type",
    rna_assay = "RNA",
    atac_assay = "ATAC",
    fragment_files = FALSE,
    metacell_args = list(),
    progress = getOption("RegCompassR.progress", TRUE)) {
  monitor <- .rc_step_monitor_start("metacells", outdir, progress)
  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)
  if (identical(fragment_files, FALSE) || is.null(fragment_files)) {
    object <- .rc_clear_signac_fragments(object, atac_assay = atac_assay)
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  pooled <- .rc_make_condition_pooled_metacells(
    object = object,
    outdir = outdir,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    fragment_files = fragment_files,
    metacell_args = metacell_args
  )
  metacell_object <- .rc_normalize_condition_metacell_object(
    pooled, rna_assay, atac_assay
  )
  if (!setequal(
    colnames(metacell_object),
    as.character(pooled$metacell_meta$metacell_id)
  )) {
    stop(
      "Merged metacell object and metadata contain different units.",
      call. = FALSE
    )
  }
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
      input_sample_col = sample_col,
      sample_col = pooled$analysis_sample_col,
      condition_col = condition_col,
      celltype_col = celltype_col,
      rna_assay = rna_assay,
      atac_assay = atac_assay,
      fragment_files = fragment_files,
      metacell_args = modifyList(list(gamma = 30L), metacell_args)
    )
  )
  class(answer) <- c("regcompass_metacell_step", "list")
  answer <- .rc_step_monitor_finish(answer, monitor)
  saveRDS(answer, file.path(outdir, "step_metacells.rds"))
  answer
}
