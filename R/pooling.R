#' Drop cells with NA grouping values
#'
#' @param meta Cell metadata.
#' @param group_cols Grouping columns that define metacell strata.
#' @return Metadata with rows containing NA in grouping columns removed.
#' @export
rc_drop_na_grouping <- function(meta, group_cols) {
  bad <- rowSums(is.na(meta[, group_cols, drop = FALSE])) > 0
  if (any(bad)) warning(sum(bad), " cells removed due to NA in grouping columns", call. = FALSE)
  meta[!bad, , drop = FALSE]
}
