# Internal runtime RDS storage policy.
#
# Package functions call `saveRDS()` without namespace qualification. Defining
# the function in the package namespace centralizes two storage rules:
#
# 1. runtime RDS files use explicit gzip compression while retaining the
#    conventional `.rds` extension;
# 2. the historical `single_cell_grn.rds` intermediate is not written because
#    the complete Stage 1 checkpoint is already stored as `step_grn.rds`.
saveRDS <- function(
    object, file = "", ascii = FALSE, version = NULL,
    compress = TRUE, refhook = NULL) {
  is_path <- is.character(file) && length(file) == 1L &&
    !is.na(file) && nzchar(file)
  if (is_path && identical(basename(file), "single_cell_grn.rds")) {
    return(invisible(NULL))
  }
  resolved_compression <- if (is_path) "gzip" else compress
  base::saveRDS(
    object = object,
    file = file,
    ascii = ascii,
    version = version,
    compress = resolved_compression,
    refhook = refhook
  )
}
