#' Validate an annotated Seurat v4 multiome object
#'
#' `rc_validate_seurat()` checks that a pre-annotated Seurat v4 object contains
#' the assays and cell metadata required by RegCompassR v0.1. It does not
#' perform clustering, WNN construction, or any other preprocessing.
#'
#' @param object A Seurat v4 object with paired RNA and ATAC assays.
#' @param rna_assay Name of the RNA assay, usually `"RNA"` or `"SCT"`.
#' @param atac_assay Name of the ATAC/chromatin assay, usually `"ATAC"`.
#' @param sample_col Metadata column containing sample identifiers.
#' @param celltype_col Metadata column containing cell type annotations.
#' @param condition_col Optional metadata column containing biological condition.
#' @param batch_col Optional metadata column containing batch labels.
#' @param embedding Optional dimensional reduction name to require.
#'
#' @return Invisibly returns `TRUE` when all requested inputs are present.
#' @export
rc_validate_seurat <- function(object,
                               rna_assay = "RNA",
                               atac_assay = "ATAC",
                               sample_col = "sample_id",
                               celltype_col = "cell_type",
                               condition_col = NULL,
                               batch_col = NULL,
                               state_col = NULL,
                               embedding = NULL) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from class 'Seurat'.", call. = FALSE)
  }

  meta <- object@meta.data
  required_meta <- c(sample_col, celltype_col, condition_col, batch_col, state_col)
  required_meta <- required_meta[!is.na(required_meta) & nzchar(required_meta)]
  missing_meta <- setdiff(required_meta, colnames(meta))
  if (length(missing_meta) > 0) {
    stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "), call. = FALSE)
  }

  assay_names <- names(object@assays)
  if (!rna_assay %in% assay_names) {
    stop("RNA assay not found: ", rna_assay, call. = FALSE)
  }
  if (!atac_assay %in% assay_names) {
    stop("ATAC assay not found: ", atac_assay, call. = FALSE)
  }

  rna_cells <- colnames(SeuratObject::GetAssayData(object = object, assay = rna_assay, slot = "counts"))
  atac_cells <- colnames(SeuratObject::GetAssayData(object = object, assay = atac_assay, slot = "counts"))
  if (!setequal(rna_cells, atac_cells)) {
    stop("RNA and ATAC assays must contain the same cell barcodes.", call. = FALSE)
  }

  if (!is.null(embedding) && !embedding %in% names(object@reductions)) {
    stop("Embedding/reduction not found: ", embedding, call. = FALSE)
  }

  invisible(TRUE)
}

#' Extract RegCompassR v0.1 inputs from an annotated Seurat v4 multiome object
#'
#' `rc_extract_inputs()` validates the requested object components and returns
#' RNA counts, ATAC counts, cell metadata, and an optional embedding in a plain
#' list for downstream RegCompassR stages.
#'
#' @inheritParams rc_validate_seurat
#' @param rna_slot Assay data slot to extract from the RNA assay. Seurat v4 uses
#' slots such as `"counts"`, `"data"`, and `"scale.data"`.
#' @param atac_slot Assay data slot to extract from the ATAC assay. Seurat v4
#' uses slots such as `"counts"` and `"data"`.
#'
#' @return A list with elements `rna`, `atac`, `meta`, and `embedding`.
#' @export
rc_extract_inputs <- function(object,
                              rna_assay = "RNA",
                              rna_slot = "counts",
                              atac_assay = "ATAC",
                              atac_slot = "counts",
                              sample_col = "sample_id",
                              celltype_col = "cell_type",
                              condition_col = NULL,
                              batch_col = NULL,
                              state_col = NULL,
                              embedding = NULL) {
  if (!identical(rna_slot, "counts") || !identical(atac_slot, "counts")) stop("Main RegCompassR workflow requires counts slots for RNA and ATAC extraction.", call. = FALSE)

  rc_validate_seurat(
    object = object,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    sample_col = sample_col,
    celltype_col = celltype_col,
    condition_col = condition_col,
    batch_col = batch_col,
    state_col = state_col,
    embedding = embedding
  )

  rna <- SeuratObject::GetAssayData(object = object, assay = rna_assay, slot = rna_slot)
  atac <- SeuratObject::GetAssayData(object = object, assay = atac_assay, slot = atac_slot)
  meta <- object@meta.data

  emb <- NULL
  if (!is.null(embedding)) {
    emb <- object@reductions[[embedding]]@cell.embeddings
  }

  list(rna_counts = rna, atac_counts = atac, rna = rna, atac = atac, meta = meta, embedding = emb)
}

#' Seurat v4 validation alias following the development plan naming
#' @export
rc_validate_seurat_v4 <- rc_validate_seurat

#' Seurat v4 extraction alias following the development plan naming
#' @export
rc_extract_seurat_v4 <- rc_extract_inputs

#' Get assay counts from a Seurat v4 object
#' @export
rc_get_assay_counts <- function(object, assay) {
  SeuratObject::GetAssayData(object = object, assay = assay, slot = "counts")
}

#' Seurat v4 validation alias following the development plan naming
#' @export
rc_validate_seurat_v4 <- rc_validate_seurat

#' Seurat v4 extraction alias following the development plan naming
#' @export
rc_extract_seurat_v4 <- rc_extract_inputs

#' Get assay counts from a Seurat v4 object
#' @export
rc_get_assay_counts <- function(object, assay) {
  SeuratObject::GetAssayData(object = object, assay = assay, slot = "counts")
}
