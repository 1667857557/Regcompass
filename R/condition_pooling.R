.rc_condition_only_sample_col <- function(meta) {
  candidate <- ".rc_condition_pool_id"
  while (candidate %in% colnames(meta)) candidate <- paste0(candidate, "_")
  candidate
}

.rc_make_condition_pooled_metacells <- function(
    object, outdir,
    sample_col = NULL,
    condition_col = "condition",
    celltype_col = "cell_type",
    rna_assay = "RNA",
    atac_assay = "ATAC",
    fragment_files = FALSE,
    metacell_args = list()) {
  if (!inherits(object, "Seurat")) stop("`object` must inherit from Seurat.", call. = FALSE)
  if (!is.list(metacell_args)) stop("`metacell_args` must be a list.", call. = FALSE)
  required <- c(condition_col, celltype_col)
  missing <- setdiff(required, colnames(object@meta.data))
  if (length(missing)) stop("Missing metadata columns: ", paste(missing, collapse = ", "), call. = FALSE)
  invalid <- vapply(object@meta.data[, required, drop = FALSE], function(x) anyNA(x) || any(!nzchar(trimws(as.character(x)))), logical(1))
  if (any(invalid)) stop("Condition and cell-type metadata must be complete.", call. = FALSE)
  if (!identical(fragment_files, FALSE) && !is.null(fragment_files)) {
    stop("The canonical condition-pooled path requires `fragment_files = FALSE` and aggregates the existing ATAC peak-count assay.", call. = FALSE)
  }
  unsupported <- intersect(names(metacell_args), c("sample_balance", "sample_balance_seed"))
  if (length(unsupported)) stop("Sample balancing is not part of the canonical workflow: ", paste(unsupported, collapse = ", "), call. = FALSE)
  if (is.null(metacell_args$gamma)) metacell_args$gamma <- 75L
  internal_sample_col <- .rc_condition_only_sample_col(object@meta.data)
  object@meta.data[[internal_sample_col]] <- paste0(as.character(object@meta.data[[condition_col]]), "__condition_pool")
  reserved <- intersect(names(metacell_args), c(
    "object", "outdir", "sample_col", "condition_col", "celltype_col",
    "rna_assay", "atac_assay", "fragment_files", "save_metacell_object",
    "save_counts", "save_fragments", "require_fragment_aggregation",
    "fragment_aggregation_backend", "on_stratum_error"
  ))
  if (length(reserved)) stop("`metacell_args` cannot override workflow fields: ", paste(reserved, collapse = ", "), call. = FALSE)
  defaults <- list(
    object = object,
    outdir = outdir,
    sample_col = internal_sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    fragment_files = FALSE,
    save_metacell_object = TRUE,
    save_counts = TRUE,
    save_fragments = FALSE,
    require_fragment_aggregation = FALSE,
    fragment_aggregation_backend = "none",
    on_stratum_error = "stop"
  )
  defaults[names(metacell_args)] <- NULL
  pooled <- do.call(rc_make_supercell2_metacells, c(defaults, metacell_args))
  meta <- pooled$metacell_meta
  if (!is.data.frame(meta) || !nrow(meta)) stop("Condition-pooled SuperCell2 produced no metacells.", call. = FALSE)
  meta$pooled_sample_id <- meta[[internal_sample_col]]
  meta$pooling_scope <- "condition_x_celltype"
  meta$sample_weighting <- "none"
  meta$sample_col_role <- "internal_condition_pool_id"
  pooled$metacell_meta <- meta
  pooled$input_sample_col <- sample_col
  pooled$analysis_sample_col <- internal_sample_col
  pooled$condition_col <- condition_col
  pooled$celltype_col <- celltype_col
  pooled$pooling_scope <- "condition_x_celltype"
  pooled$sample_weighting <- "none"
  pooled$input_design <- list(
    metacell_grouping = c(condition_col, celltype_col),
    condition_only_stratification = TRUE,
    gamma = metacell_args$gamma,
    inference_policy = "cells are stratified by condition and cell type; sample metadata are not used for selection, weighting or grouping"
  )
  pooled
}
