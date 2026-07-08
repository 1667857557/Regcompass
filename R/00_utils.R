# Utility helpers shared across RegCompassR modules.
NULL

.rc_first_existing_col <- function(candidates, x, fallback = NULL) {
  hit <- intersect(candidates, colnames(x))
  if (length(hit) > 0L && !is.na(hit[[1]]) && nzchar(hit[[1]])) return(hit[[1]])
  if (!is.null(fallback)) return(fallback)
  colnames(x)[1]
}

#' Drop rows with missing grouping values
#'
#' Shared metadata cleanup for workflows that split objects by sample, condition,
#' cell type, or other grouping columns.
#' @export
rc_drop_na_grouping <- function(meta, group_cols) {
  if (!is.data.frame(meta)) stop("`meta` must be a data.frame.", call. = FALSE)
  group_cols <- group_cols[!is.na(group_cols) & nzchar(group_cols)]
  missing_cols <- setdiff(group_cols, colnames(meta))
  if (length(missing_cols) > 0L) stop("`meta` is missing grouping columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  if (length(group_cols) == 0L) return(meta)
  keep <- stats::complete.cases(meta[, group_cols, drop = FALSE])
  meta[keep, , drop = FALSE]
}
