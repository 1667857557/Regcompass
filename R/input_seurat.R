.rc_get_assay_counts <- function(object, assay) {
  seurat_object_version <- tryCatch(
    utils::packageVersion("SeuratObject"),
    error = function(e) utils::package_version("0.0.0")
  )
  if (seurat_object_version >= utils::package_version("5.0.0")) {
    SeuratObject::GetAssayData(
      object = object,
      assay = assay,
      layer = "counts"
    )
  } else {
    SeuratObject::GetAssayData(
      object = object,
      assay = assay,
      slot = "counts"
    )
  }
}
