.rc_celltype_col_from_pooled <- function(pooled) {
  design <- pooled$input_design$condition_celltype_sample_count %||% NULL
  if (is.data.frame(design) && ncol(design) >= 2L) {
    candidate <- colnames(design)[[2L]]
    if (!is.na(candidate) && nzchar(candidate)) return(candidate)
  }
  if ("cell_type" %in% colnames(pooled$metacell_meta)) return("cell_type")
  stop(
    "Unable to infer the cell-type column for shared ATAC TF-IDF normalization.",
    call. = FALSE
  )
}

.rc_drop_zero_count_atac_features <- function(
    object, atac_assay = "ATAC", context = "ATAC normalization") {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  counts <- .rc_get_assay_counts(object, atac_assay)
  peak_totals <- Matrix::rowSums(counts)
  keep <- is.finite(peak_totals) & peak_totals > 0
  if (!any(keep)) {
    stop(context, " has no peaks with a positive total count.", call. = FALSE)
  }
  diagnostics <- list(
    n_input_peaks = nrow(counts),
    n_zero_count_peaks_excluded = sum(!keep),
    n_retained_peaks = sum(keep),
    zero_count_peak_policy = "exclude_before_tfidf_and_pando"
  )
  if (any(!keep)) {
    assay_object <- object[[atac_assay]]
    assay_object <- subset(
      assay_object,
      features = rownames(counts)[keep]
    )
    object[[atac_assay]] <- assay_object
  }
  list(object = object, diagnostics = diagnostics)
}

.rc_apply_celltype_shared_tfidf <- function(
    object, celltype_col, atac_assay = "ATAC",
    method = 1, scale.factor = 1e4) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  if (!celltype_col %in% colnames(object@meta.data)) {
    stop("Missing cell-type metadata column: ", celltype_col, call. = FALSE)
  }
  filtered <- .rc_drop_zero_count_atac_features(
    object,
    atac_assay = atac_assay,
    context = "Cell-type-shared TF-IDF"
  )
  object <- filtered$object
  counts <- .rc_get_assay_counts(object, atac_assay)
  units <- colnames(counts)
  meta_index <- match(units, rownames(object@meta.data))
  if (anyNA(meta_index)) {
    stop("ATAC metacells do not align with object metadata.", call. = FALSE)
  }
  meta <- object@meta.data[meta_index, , drop = FALSE]
  cell_type <- trimws(as.character(meta[[celltype_col]]))
  if (anyNA(cell_type) || any(!nzchar(cell_type))) {
    stop("Cell-type metadata are incomplete for shared TF-IDF.", call. = FALSE)
  }

  groups <- split(units, cell_type)
  normalized <- lapply(groups, function(group_units) {
    Signac::RunTFIDF(
      object = counts[, group_units, drop = FALSE],
      method = method,
      scale.factor = scale.factor,
      verbose = FALSE
    )
  })
  tfidf <- do.call(cbind, normalized)
  tfidf <- tfidf[rownames(counts), units, drop = FALSE]
  if (!identical(dim(tfidf), dim(counts)) ||
      !identical(dimnames(tfidf), dimnames(counts))) {
    stop("Cell-type-shared TF-IDF did not preserve the ATAC matrix layout.",
         call. = FALSE)
  }

  assay_object <- object[[atac_assay]]
  assay_object <- SeuratObject::SetAssayData(
    object = assay_object,
    slot = "data",
    new.data = tfidf
  )
  object[[atac_assay]] <- assay_object
  object@misc$regcompass_atac_normalization <- c(
    list(
      method = "Signac_TFIDF",
      scope = "cell_type_across_conditions",
      celltype_col = celltype_col,
      idf_reference = "all condition-pooled metacells of the same cell type",
      n_metacells_by_celltype = vapply(groups, length, integer(1)),
      tfidf_method = method,
      scale_factor = scale.factor
    ),
    filtered$diagnostics
  )
  object
}

# Canonical correction: estimate one shared IDF reference per cell type across
# all conditions, then retain those normalized values for Pando and Layer 1.
.rc_normalize_condition_metacell_object_v170 <- function(
    pooled, rna_assay = "RNA", atac_assay = "ATAC") {
  object <- rc_load_or_merge_metacell_objects(
    pooled$metacell_objects,
    fragment_manifest = NULL,
    metacell_meta = pooled$metacell_meta,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    require_complete_fragments = FALSE
  )
  object <- Seurat::NormalizeData(object, assay = rna_assay, verbose = FALSE)
  object <- .rc_apply_celltype_shared_tfidf(
    object,
    celltype_col = .rc_celltype_col_from_pooled(pooled),
    atac_assay = atac_assay
  )
  object
}
.rc_normalize_condition_metacell_object <-
  .rc_normalize_condition_metacell_object_v170
