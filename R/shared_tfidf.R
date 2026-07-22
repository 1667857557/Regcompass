.rc_drop_zero_count_atac_features <- function(object, atac_assay = "ATAC", context = "ATAC normalization") {
  if (!inherits(object, "Seurat")) stop("`object` must inherit from Seurat.", call. = FALSE)
  counts <- .rc_get_assay_counts(object, atac_assay)
  totals <- Matrix::rowSums(counts)
  keep <- is.finite(totals) & totals > 0
  if (!any(keep)) stop(context, " has no peaks with a positive total count.", call. = FALSE)
  diagnostics <- list(
    n_input_peaks = nrow(counts),
    n_zero_count_peaks_excluded = sum(!keep),
    n_retained_peaks = sum(keep),
    zero_count_peak_policy = "exclude_before_tfidf_and_pando"
  )
  if (any(!keep)) object[[atac_assay]] <- subset(object[[atac_assay]], features = rownames(counts)[keep])
  list(object = object, diagnostics = diagnostics)
}

.rc_apply_celltype_shared_tfidf <- function(object, celltype_col, atac_assay = "ATAC", method = 1, scale.factor = 1e4) {
  if (!inherits(object, "Seurat")) stop("`object` must inherit from Seurat.", call. = FALSE)
  if (!celltype_col %in% colnames(object@meta.data)) stop("Missing cell-type metadata column: ", celltype_col, call. = FALSE)
  filtered <- .rc_drop_zero_count_atac_features(object, atac_assay, "Cell-type-shared TF-IDF")
  object <- filtered$object
  counts <- .rc_get_assay_counts(object, atac_assay)
  units <- colnames(counts)
  meta <- object@meta.data[match(units, rownames(object@meta.data)), , drop = FALSE]
  cell_type <- trimws(as.character(meta[[celltype_col]]))
  if (anyNA(cell_type) || any(!nzchar(cell_type))) stop("Cell-type metadata are incomplete for shared TF-IDF.", call. = FALSE)
  groups <- split(units, cell_type)
  normalized <- lapply(groups, function(group_units) {
    Signac::RunTFIDF(counts[, group_units, drop = FALSE], method = method, scale.factor = scale.factor, verbose = FALSE)
  })
  tfidf <- do.call(cbind, normalized)
  tfidf <- tfidf[rownames(counts), units, drop = FALSE]
  if (!identical(dimnames(tfidf), dimnames(counts))) stop("Cell-type-shared TF-IDF did not preserve the ATAC matrix layout.", call. = FALSE)
  assay_object <- object[[atac_assay]]
  assay_object <- SeuratObject::SetAssayData(assay_object, slot = "data", new.data = tfidf)
  object[[atac_assay]] <- assay_object
  object@misc$regcompass_atac_normalization <- c(list(
    method = "Signac_TFIDF",
    scope = "cell_type_across_conditions",
    celltype_col = celltype_col,
    idf_reference = "all single cells or metacells of the same cell type across conditions",
    n_units_by_celltype = vapply(groups, length, integer(1)),
    tfidf_method = method,
    scale_factor = scale.factor
  ), filtered$diagnostics)
  object
}

.rc_normalize_single_cell_grn_object <- function(object, condition_col = "condition", celltype_col = "cell_type", rna_assay = "RNA", atac_assay = "ATAC") {
  if (!inherits(object, "Seurat")) stop("`object` must inherit from Seurat.", call. = FALSE)
  required <- c(condition_col, celltype_col)
  missing <- setdiff(required, colnames(object@meta.data))
  if (length(missing)) stop("Missing metadata columns: ", paste(missing, collapse = ", "), call. = FALSE)
  invalid <- vapply(object@meta.data[, required, drop = FALSE], function(x) anyNA(x) || any(!nzchar(trimws(as.character(x)))), logical(1))
  if (any(invalid)) stop("Condition and cell-type metadata must be complete.", call. = FALSE)
  object <- Seurat::NormalizeData(object, assay = rna_assay, verbose = FALSE)
  object <- .rc_apply_celltype_shared_tfidf(object, celltype_col = celltype_col, atac_assay = atac_assay)
  object@misc$regcompass_grn_normalization <- list(
    rna = "global_single_cell_NormalizeData",
    atac = "cell_type_shared_TFIDF_across_conditions",
    condition_col = condition_col,
    celltype_col = celltype_col
  )
  object
}

.rc_normalize_condition_metacell_object <- function(pooled, rna_assay = "RNA", atac_assay = "ATAC") {
  object <- rc_load_or_merge_metacell_objects(
    pooled$metacell_objects,
    fragment_manifest = NULL,
    metacell_meta = pooled$metacell_meta,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    require_complete_fragments = FALSE
  )
  object <- Seurat::NormalizeData(object, assay = rna_assay, verbose = FALSE)
  .rc_apply_celltype_shared_tfidf(object, celltype_col = pooled$celltype_col, atac_assay = atac_assay)
}
