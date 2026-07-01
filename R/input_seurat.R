#' Validate an annotated Seurat v4 RNA+ATAC counts object
#'
#' RegCompassR Layer 1 starts from raw counts. The object must contain paired RNA
#' and ATAC assays with identical cell barcodes plus sample and cell-type labels.
#' Optional condition labels define biological strata, but clustering/WNN is not
#' rerun inside RegCompassR.
#' @export
rc_validate_seurat_v4 <- function(object,
                                  rna_assay = "RNA",
                                  atac_assay = "ATAC",
                                  sample_col = "sample_id",
                                  celltype_col = "cell_type",
                                  condition_col = NULL) {
  if (!inherits(object, "Seurat")) stop("`object` must inherit from class 'Seurat'.", call. = FALSE)

  meta <- object@meta.data
  required_meta <- c(sample_col, celltype_col, condition_col)
  required_meta <- required_meta[!is.null(required_meta) & !is.na(required_meta) & nzchar(required_meta)]
  missing_meta <- setdiff(required_meta, colnames(meta))
  if (length(missing_meta) > 0L) stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "), call. = FALSE)

  assay_names <- names(object@assays)
  if (!rna_assay %in% assay_names) stop("RNA assay not found: ", rna_assay, call. = FALSE)
  if (!atac_assay %in% assay_names) stop("ATAC assay not found: ", atac_assay, call. = FALSE)

  rna_counts <- SeuratObject::GetAssayData(object = object, assay = rna_assay, slot = "counts")
  atac_counts <- SeuratObject::GetAssayData(object = object, assay = atac_assay, slot = "counts")
  if (is.null(dim(rna_counts)) || is.null(dim(atac_counts))) stop("RNA and ATAC assays must contain count matrices.", call. = FALSE)
  if (!setequal(colnames(rna_counts), colnames(atac_counts))) stop("RNA and ATAC assays must contain the same cell barcodes.", call. = FALSE)
  if (any(Matrix::colSums(rna_counts) < 0) || any(Matrix::colSums(atac_counts) < 0)) stop("Counts must be non-negative.", call. = FALSE)

  invisible(TRUE)
}

#' Extract RNA/ATAC counts and metadata from an annotated Seurat v4 object
#' @export
rc_extract_seurat_v4 <- function(object,
                                 rna_assay = "RNA",
                                 atac_assay = "ATAC",
                                 sample_col = "sample_id",
                                 celltype_col = "cell_type",
                                 condition_col = NULL) {
  rc_validate_seurat_v4(object, rna_assay, atac_assay, sample_col, celltype_col, condition_col)
  rna_counts <- SeuratObject::GetAssayData(object = object, assay = rna_assay, slot = "counts")
  atac_counts <- SeuratObject::GetAssayData(object = object, assay = atac_assay, slot = "counts")
  list(rna_counts = rna_counts, atac_counts = atac_counts, meta = object@meta.data)
}

#' Backward-compatible aliases
#' @export
rc_validate_seurat <- rc_validate_seurat_v4

#' @export
rc_extract_inputs <- rc_extract_seurat_v4

#' Get raw counts from a Seurat v4 assay
#' @export
rc_get_assay_counts <- function(object, assay) {
  SeuratObject::GetAssayData(object = object, assay = assay, slot = "counts")
}
