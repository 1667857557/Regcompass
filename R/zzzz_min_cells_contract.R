# Final Stage 1 threshold contract. This file is collated after the workflow
# hardening overrides so one fixed `min_cells` value drives both prefiltering and
# the downstream standard/condition Pando calls.

.rc_resolve_stage1_min_cells_contract <- function(pando_args) {
  if (!is.list(pando_args)) {
    stop("`pando_args` must be a list.", call. = FALSE)
  }
  supplied <- pando_args$min_cells %||% .rc_stage1_min_cells_fixed
  supplied <- suppressWarnings(as.integer(supplied[[1L]]))
  if (!is.finite(supplied) || supplied != .rc_stage1_min_cells_fixed) {
    message("Stage 1 `min_cells` is fixed at 300; overriding the supplied value.")
  }
  pando_args$min_cells <- .rc_stage1_min_cells_fixed
  list(min_cells = pando_args$min_cells, pando_args = pando_args)
}

.rc_filter_stage1_groups_by_min_cells <- function(
    object, celltype_col, cell_type, min_cells) {
  .rc_prefilter_stage1_celltypes(
    object = object,
    celltype_col = celltype_col,
    cell_type = cell_type,
    min_cells = min_cells
  )
}

#' Infer regulatory evidence using the fixed Stage 1 `min_cells` contract
#'
#' Stage 1 fixes `pando_args$min_cells` at 300. The same value is first used to
#' remove broad cell types below the threshold before normalization and is then
#' passed unchanged into the standard or condition-aware Pando runtime.
#' Globally zero ATAC peaks are removed before TF-IDF or motif analysis. By
#' default stale Signac fragment references are cleared because Stage 1 uses the
#' in-memory peak matrix and genome sequence rather than fragment files.
#'
#' @param fragment_files Preserve existing Signac fragment references when TRUE.
#' The default FALSE clears them before Stage 1.
#' @export
rc_regcompass_step_grn <- function(
    object, gem, outdir, genome,
    pfm = NULL,
    species = c("auto", "human", "mouse"),
    condition_col = "condition",
    celltype_col = "cell_type",
    cell_type = NULL,
    rna_assay = "RNA",
    atac_assay = "ATAC",
    fragment_files = FALSE,
    pando_args = list(),
    parallel = TRUE,
    BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  threshold_contract <- .rc_resolve_stage1_min_cells_contract(pando_args)
  pando_args <- threshold_contract$pando_args
  min_cells <- threshold_contract$min_cells

  preserve_fragments <- .rc_validate_stage1_fragment_policy(fragment_files)
  filtered_groups <- .rc_filter_stage1_groups_by_min_cells(
    object = object,
    celltype_col = celltype_col,
    cell_type = cell_type,
    min_cells = min_cells
  )
  object <- filtered_groups$object
  filtered_groups$diagnostics$threshold_source <-
    "fixed_pando_args_min_cells"
  object@misc$regcompass_stage1_group_filter <- filtered_groups$diagnostics
  object@misc$regcompass_stage1_min_cells_contract <- list(
    min_cells = min_cells,
    source = "pando_args$min_cells",
    fixed = TRUE,
    applied_before_normalization = TRUE,
    passed_to_pando = TRUE
  )

  if (!preserve_fragments) {
    object <- .rc_clear_signac_fragments(object, atac_assay = atac_assay)
  }
  zero_filtered <- .rc_drop_zero_count_atac_features(
    object, atac_assay, "Stage 1 min_cells prefilter"
  )
  object <- zero_filtered$object
  object@misc$regcompass_stage1_zero_peak_filter <- zero_filtered$diagnostics
  object@misc$regcompass_stage1_fragment_policy <- list(
    fragment_files = preserve_fragments,
    policy = if (preserve_fragments) "preserve" else "clear_before_stage1"
  )

  motif_args <- pando_args$pando_motif_args %||% list()
  if (!is.list(motif_args)) {
    stop("`pando_motif_args` must be a list.", call. = FALSE)
  }
  if (.rc_pando_supports_motif_cache()) {
    motif_args$cache_dir <- motif_args$cache_dir %||%
      file.path(outdir, "motif_cache")
    motif_args$reuse_cache <- motif_args$reuse_cache %||% TRUE
  }
  pando_args$pando_motif_args <- motif_args

  duplicate_file <- file.path(outdir, "single_cell_grn.rds")
  if (file.exists(duplicate_file)) unlink(duplicate_file, force = TRUE)
  on.exit({
    if (file.exists(duplicate_file)) unlink(duplicate_file, force = TRUE)
  }, add = TRUE)

  .rc_original_step_grn_hardening(
    object = object,
    gem = gem,
    outdir = outdir,
    genome = genome,
    pfm = pfm,
    species = species,
    condition_col = condition_col,
    celltype_col = celltype_col,
    cell_type = if (is.null(cell_type)) NULL else
      filtered_groups$retained_cell_types,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    pando_args = pando_args,
    parallel = parallel,
    BPPARAM = BPPARAM,
    progress = progress
  )
}
