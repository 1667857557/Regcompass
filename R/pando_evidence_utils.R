.rc_pando_assay_data <- function(object, assay) {
  value <- tryCatch(
    .rc_get_assay_matrix(object, assay, "data"),
    error = function(e) NULL
  )
  if (is.null(value) || nrow(value) == 0L || ncol(value) == 0L) {
    stop(
      "Pando evidence projection requires one non-empty normalized `data` ",
      "layer per assay.",
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
