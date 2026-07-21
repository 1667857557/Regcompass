.rc_condition_only_sample_col <- function(meta) {
  candidate <- ".rc_condition_only_pool_id"
  while (candidate %in% colnames(meta)) candidate <- paste0(candidate, "_")
  candidate
}

.rc_prepare_condition_only_object <- function(object, condition_col) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  if (!is.character(condition_col) || length(condition_col) != 1L ||
      is.na(condition_col) || !nzchar(condition_col)) {
    stop("`condition_col` must name one metadata column.", call. = FALSE)
  }
  if (!condition_col %in% colnames(object@meta.data)) {
    stop("Missing metadata column: ", condition_col, call. = FALSE)
  }
  condition <- trimws(as.character(object@meta.data[[condition_col]]))
  if (anyNA(condition) || any(!nzchar(condition))) {
    stop("Condition metadata must be complete and non-empty.", call. = FALSE)
  }
  internal_sample_col <- .rc_condition_only_sample_col(object@meta.data)
  object@meta.data[[internal_sample_col]] <- paste0(
    condition,
    "__condition_pool"
  )
  list(object = object, sample_col = internal_sample_col)
}

.rc_make_condition_pooled_metacells_unbalanced <-
  .rc_make_condition_pooled_metacells

.rc_make_condition_pooled_metacells_v170 <- function(
    object, outdir,
    sample_col = NULL,
    condition_col = "condition",
    celltype_col = "cell_type",
    rna_assay = "RNA",
    atac_assay = "ATAC",
    fragment_files = FALSE,
    metacell_args = list(),
    strict_biological_defaults = TRUE) {
  if (!is.list(metacell_args)) {
    stop("`metacell_args` must be a list.", call. = FALSE)
  }
  if (!is.null(sample_col) &&
      (!is.character(sample_col) || length(sample_col) != 1L ||
       is.na(sample_col) || !nzchar(sample_col))) {
    stop("`sample_col` must be NULL or one metadata-column name.", call. = FALSE)
  }

  ignored_balance_args <- intersect(
    names(metacell_args),
    c("sample_balance", "sample_balance_seed")
  )
  if (length(ignored_balance_args)) {
    warning(
      paste0(
        "Ignoring `", paste(ignored_balance_args, collapse = "`, `"),
        "`: metacells are stratified only by condition and cell type; ",
        "biological-sample labels do not alter cell selection, weighting, or grouping."
      ),
      call. = FALSE
    )
    metacell_args[ignored_balance_args] <- NULL
  }

  prepared <- .rc_prepare_condition_only_object(object, condition_col)
  pooled <- withCallingHandlers(
    .rc_make_condition_pooled_metacells_unbalanced(
      object = prepared$object,
      outdir = outdir,
      sample_col = prepared$sample_col,
      condition_col = condition_col,
      celltype_col = celltype_col,
      rna_assay = rna_assay,
      atac_assay = atac_assay,
      fragment_files = fragment_files,
      metacell_args = metacell_args,
      strict_biological_defaults = FALSE
    ),
    warning = function(warning) {
      if (grepl(
        "Condition-pooled analysis requires at least two biological samples",
        conditionMessage(warning),
        fixed = TRUE
      )) {
        invokeRestart("muffleWarning")
      }
    }
  )

  sample_summary_cols <- intersect(
    c(
      "n_biological_samples", "dominant_sample_fraction",
      "effective_sample_n"
    ),
    colnames(pooled$metacell_meta)
  )
  for (column in sample_summary_cols) {
    pooled$metacell_meta[[column]] <- if (column == "n_biological_samples") {
      NA_integer_
    } else {
      NA_real_
    }
  }
  if ("samples_mixed_within_condition" %in% colnames(pooled$metacell_meta)) {
    pooled$metacell_meta$samples_mixed_within_condition <- FALSE
  }

  pooled$metacell_meta$sample_weighting <-
    "not_applicable_condition_only_stratification"
  pooled$metacell_meta$sample_balance <- FALSE
  pooled$metacell_meta$sample_col_role <-
    "internal_condition_pool_id_not_biological_sample"
  if (is.data.frame(pooled$sample_composition)) {
    pooled$sample_composition$sample_col_role <-
      "internal_condition_pool_id_not_biological_sample"
  }
  if (is.data.frame(pooled$sample_composition_summary)) {
    pooled$sample_composition_summary$sample_provenance_available <- FALSE
  }

  pooled$input_sample_col <- sample_col
  pooled$analysis_sample_col <- prepared$sample_col
  pooled$sample_balance <- FALSE
  pooled$sample_balance_seed <- NULL
  pooled$sample_weighting <- "not_applicable_condition_only_stratification"
  pooled$sample_balance_diagnostics <- data.frame()
  pooled$sample_balance_summary <- list(
    sample_balance = FALSE,
    sample_balance_seed = NULL,
    sample_weighting = pooled$sample_weighting,
    n_input_cells = ncol(object),
    n_retained_cells = ncol(object),
    n_excluded_cells = 0L
  )
  pooled$pooling_scope <- "condition_x_celltype"
  pooled$input_design$sample_col_role <-
    "ignored_not_used_for_stratification_weighting_or_cell_selection"
  pooled$input_design$metacell_grouping <- c(condition_col, celltype_col)
  pooled$input_design$condition_only_stratification <- TRUE
  pooled$input_design$sample_balance <- pooled$sample_balance_summary
  pooled$input_design$inference_policy <- paste(
    "cells are stratified only by condition and cell type before metacell",
    "construction; sample metadata are not used"
  )
  pooled
}

.rc_make_condition_pooled_metacells <-
  .rc_make_condition_pooled_metacells_v170
