.rc_clear_signac_fragments <- function(object, atac_assay = "ATAC") {
  if (!inherits(object, "Seurat")) return(object)
  if (!requireNamespace("Signac", quietly = TRUE)) return(object)
  if (!atac_assay %in% names(object@assays)) return(object)
  if (!inherits(object[[atac_assay]], "ChromatinAssay")) return(object)
  fragment_setter <- get("Fragments<-", envir = asNamespace("Signac"))
  object[[atac_assay]] <- fragment_setter(
    object[[atac_assay]], value = list()
  )
  object
}
