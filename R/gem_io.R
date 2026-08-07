# Internal GEM persistence helpers ---------------------------------------------

.rc_validate_model_info <- function(model_info) {
  if (is.null(model_info) || !is.list(model_info)) {
    stop("GEM is missing `model_info`.", call. = FALSE)
  }
  required <- c(
    "source", "version", "commit", "checksum", "conversion_date"
  )
  present <- vapply(required, function(name) {
    !is.null(model_info[[name]]) &&
      length(model_info[[name]]) > 0L &&
      !all(is.na(model_info[[name]]))
  }, logical(1))
  if (!all(present)) {
    stop(
      "`model_info` is incomplete: ",
      paste(required[!present], collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.rc_read_gem <- function(file) {
  gem <- readRDS(file)
  rc_validate_gem(gem)
  .rc_validate_model_info(gem$model_info)
  gem
}
