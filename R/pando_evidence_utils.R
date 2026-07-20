.rc_pando_assay_data <- function(object, assay) {
  value <- tryCatch(
    SeuratObject::GetAssayData(object, assay = assay, slot = "data"),
    error = function(e) NULL
  )
  if (is.null(value)) {
    value <- tryCatch(
      SeuratObject::GetAssayData(object, assay = assay, layer = "data"),
      error = function(e) NULL
    )
  }
  if (is.null(value) || nrow(value) == 0L || ncol(value) == 0L) {
    stop(
      "Pando evidence projection requires non-empty normalized assay data.",
      call. = FALSE
    )
  }
  value
}

.rc_pando_region_key <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^([^:]+):(\\d+)-(\\d+)$", "\\1-\\2-\\3", x)
  x <- sub("^([^:]+):(\\d+):(\\d+)$", "\\1-\\2-\\3", x)
  x
}
