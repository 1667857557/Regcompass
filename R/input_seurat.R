.rc_get_assay_counts <- function(object, assay) {
  get_assay_data <- SeuratObject::GetAssayData
  args <- list(object = object, assay = assay)
  if ("layer" %in% names(formals(get_assay_data))) {
    args$layer <- "counts"
  } else {
    args$slot <- "counts"
  }
  do.call(get_assay_data, args)
}
