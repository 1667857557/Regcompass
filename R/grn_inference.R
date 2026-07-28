.rc_require_normalized_assay <- function(object, assay, context) {
  value <- .rc_pando_assay_data(object, assay)
  expected <- as.character(colnames(object))
  observed <- as.character(colnames(value))
  if (anyNA(expected) || any(!nzchar(expected)) || anyDuplicated(expected) ||
      anyNA(observed) || any(!nzchar(observed)) || anyDuplicated(observed)) {
    stop(
      context,
      " normalized assay cell identifiers must be unique and non-empty.",
      call. = FALSE
    )
  }
  missing <- setdiff(expected, observed)
  extra <- setdiff(observed, expected)
  if (length(missing) || length(extra)) {
    stop(
      context,
      " normalized assay contains different cells from the analysis object. ",
      "Missing: ", paste(utils::head(missing, 10L), collapse = ", "),
      "; extra: ", paste(utils::head(extra, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}
