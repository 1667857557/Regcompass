# Internal utility helpers loaded before the remaining package files.

.rc_first_existing_col <- function(x, candidates, fallback = NULL) {
  if (is.null(colnames(x))) {
    if (!is.null(fallback)) return(fallback)
    return(NULL)
  }
  candidates <- unique(as.character(candidates))
  hit <- intersect(candidates, colnames(x))
  if (length(hit)) return(hit[[1L]])
  if (!is.null(fallback) && length(fallback) == 1L &&
      !is.na(fallback) && nzchar(fallback)) {
    return(fallback)
  }
  NULL
}

#' Drop rows with missing grouping metadata
#'
#' @param meta Data frame containing metadata rows.
#' @param grouping_cols Character vector of grouping columns.
#' @return Filtered metadata with an attribute containing per-column drop counts.
#' @export
rc_drop_na_grouping <- function(meta, grouping_cols) {
  missing_cols <- setdiff(grouping_cols, colnames(meta))
  if (length(missing_cols) > 0L) {
    stop("Missing grouping columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  drop_by_col <- vapply(grouping_cols, function(column) {
    values <- meta[[column]]
    sum(is.na(values) | !nzchar(trimws(as.character(values))))
  }, integer(1))
  keep <- stats::complete.cases(meta[, grouping_cols, drop = FALSE])
  for (column in grouping_cols) {
    keep <- keep & nzchar(trimws(as.character(meta[[column]])))
  }
  out <- meta[keep, , drop = FALSE]
  attr(out, "dropped_na_by_column") <- drop_by_col
  out
}
