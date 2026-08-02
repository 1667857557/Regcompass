# Final cross-stage cell-set contract. Loaded after the Stage 1 min-cells
# override so Stage 1 and Stage 2 consume the same cells.

.rc_original_step_grn_cell_set_contract <- rc_regcompass_step_grn
.rc_original_step_metacells_cell_set_contract <- rc_regcompass_step_metacells

.rc_build_stage_analysis_cell_set <- function(
    object, condition_col, celltype_col, cell_type, pando_args) {
  threshold <- .rc_resolve_stage1_min_cells_contract(pando_args)
  groups <- .rc_filter_stage1_groups_by_min_cells(
    object = object,
    condition_col = condition_col,
    celltype_col = celltype_col,
    cell_type = cell_type,
    min_cells = threshold$min_cells
  )
  analysis_types <- unique(as.character(groups$pando_cell_types))
  observed_type <- trimws(as.character(
    groups$object@meta.data[[celltype_col]]
  ))
  keep <- observed_type %in% analysis_types
  analysis_cells <- rownames(groups$object@meta.data)[keep]
  if (!length(analysis_cells)) {
    stop("No cells remain for Stage 1 Pando analysis.", call. = FALSE)
  }
  analysis_object <- subset(groups$object, cells = analysis_cells)
  list(
    object = analysis_object,
    retained_cells = analysis_cells,
    retained_cell_types = analysis_types,
    diagnostics = groups$diagnostics,
    analysis_mode = groups$analysis_mode,
    condition_levels = groups$condition_levels,
    min_cells = threshold$min_cells,
    pando_args = threshold$pando_args,
    skipped_condition_cell_types = groups$skipped_condition_cell_types
  )
}

.rc_validate_stage1_cell_set <- function(grn) {
  if (!inherits(grn, "regcompass_grn_step")) {
    stop("`grn` must be the output of `rc_regcompass_step_grn()`.",
         call. = FALSE)
  }
  contract <- grn$cell_filter
  if (!is.list(contract)) {
    stop(
      "The Stage 1 result has no cell-set contract. Rerun Stage 1 with the current RegCompassR version.",
      call. = FALSE
    )
  }
  cells <- as.character(contract$retained_cells)
  types <- as.character(contract$retained_cell_types)
  if (!length(cells) || anyNA(cells) || any(!nzchar(cells)) ||
      anyDuplicated(cells)) {
    stop("The Stage 1 retained-cell contract is invalid.", call. = FALSE)
  }
  if (!length(types) || anyNA(types) || any(!nzchar(types))) {
    stop("The Stage 1 retained-cell-type contract is invalid.", call. = FALSE)
  }
  contract$retained_cells <- cells
  contract$retained_cell_types <- unique(types)
  contract
}

.rc_subset_to_stage1_cell_set <- function(object, contract) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  available <- colnames(object)
  missing <- setdiff(contract$retained_cells, available)
  if (length(missing)) {
    stop(
      "Stage 2 input is missing ", length(missing),
      " cell(s) retained by Stage 1; first missing ID: ", missing[[1L]],
      call. = FALSE
    )
  }
  filtered <- subset(object, cells = contract$retained_cells)
  if (!identical(colnames(filtered), contract$retained_cells)) {
    filtered <- filtered[, contract$retained_cells]
  }
  if (!identical(colnames(filtered), contract$retained_cells)) {
    stop("Stage 2 could not reproduce the ordered Stage 1 cell set.",
         call. = FALSE)
  }
  filtered@misc$regcompass_cross_stage_cell_set <- list(
    source = contract$source %||% "stage1_grn_result",
    n_cells = length(contract$retained_cells),
    retained_cell_types = contract$retained_cell_types,
    min_cells = contract$min_cells %||% .rc_stage1_min_cells_fixed
  )
  filtered
}

# Infer regulatory evidence and persist the exact downstream cell set.
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
  n_input <- ncol(object)
  cell_set <- .rc_build_stage_analysis_cell_set(
    object = object,
    condition_col = condition_col,
    celltype_col = celltype_col,
    cell_type = cell_type,
    pando_args = pando_args
  )
  answer <- .rc_original_step_grn_cell_set_contract(
    object = cell_set$object,
    gem = gem,
    outdir = outdir,
    genome = genome,
    pfm = pfm,
    species = species,
    condition_col = condition_col,
    celltype_col = celltype_col,
    cell_type = cell_set$retained_cell_types,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    fragment_files = fragment_files,
    pando_args = cell_set$pando_args,
    parallel = parallel,
    BPPARAM = BPPARAM,
    progress = progress
  )
  answer$cell_filter <- list(
    source = "stage1_min_cells_before_normalization",
    min_cells = cell_set$min_cells,
    n_input_cells = as.integer(n_input),
    n_retained_cells = as.integer(length(cell_set$retained_cells)),
    n_removed_cells = as.integer(n_input - length(cell_set$retained_cells)),
    retained_cells = cell_set$retained_cells,
    retained_cell_types = cell_set$retained_cell_types,
    skipped_condition_cell_types = cell_set$skipped_condition_cell_types,
    diagnostics = cell_set$diagnostics,
    analysis_mode = cell_set$analysis_mode,
    condition_levels = cell_set$condition_levels
  )
  answer$params$cell_type <- cell_set$retained_cell_types
  answer$params$n_input_cells <- as.integer(n_input)
  answer$params$n_stage_cells <- as.integer(length(cell_set$retained_cells))
  answer$params$cell_set_contract <- "stage1_exact_cell_ids_v1"
  saveRDS(answer, file.path(outdir, "step_grn.rds"))
  answer
}

# Build metacells from the exact Stage 1 cell set. `grn` is optional for
# backward compatibility; supplying it enables exact ID and parameter checks.
rc_regcompass_step_metacells <- function(
    object, outdir,
    condition_col = "condition",
    celltype_col = "cell_type",
    cell_type = NULL,
    rna_assay = "RNA",
    atac_assay = "ATAC",
    fragment_files = FALSE,
    metacell_args = list(),
    progress = getOption("RegCompassR.progress", TRUE),
    grn = NULL) {
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
      skipped_condition_cell_types = cell_set$skipped_condition_cell_types,
      diagnostics = cell_set$diagnostics,
      analysis_mode = cell_set$analysis_mode,
      condition_levels = cell_set$condition_levels
    )
  } else {
    contract <- .rc_validate_stage1_cell_set(grn)
    expected_condition_col <-
      grn$params$requested_condition_col %||% grn$params$condition_col
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
  answer <- .rc_original_step_metacells_cell_set_contract(
    object = object,
    outdir = outdir,
    condition_col = condition_col,
    celltype_col = celltype_col,
    cell_type = contract$retained_cell_types,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    fragment_files = fragment_files,
    metacell_args = metacell_args,
    progress = progress
  )
  answer$cell_filter <- contract
  answer$cell_filter$stage2_n_input_cells <- as.integer(n_input)
  answer$cell_filter$n_stage_cells <-
    as.integer(length(contract$retained_cells))
  answer$cell_filter$exact_stage1_match <- !is.null(grn)
  answer$params$n_input_cells <- as.integer(n_input)
  answer$params$n_stage_cells <- as.integer(length(contract$retained_cells))
  answer$params$cell_set_contract <- if (is.null(grn)) {
    "independent_stage1_filter_reapplication_v1"
  } else {
    "stage1_exact_cell_ids_v1"
  }
  saveRDS(answer, file.path(outdir, "step_metacells.rds"))
  answer
}
