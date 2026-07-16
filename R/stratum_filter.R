.rc_strict_stratum_cols <- function(sample_col, condition_col, celltype_col) {
  cols <- unique(c(condition_col, sample_col, celltype_col))
  cols <- cols[!is.na(cols) & nzchar(cols)]
  if (length(cols) != 3L) {
    stop("Condition, sample and cell-type columns must be distinct.", call. = FALSE)
  }
  cols
}

.rc_add_stratum_id <- function(meta, cols) {
  missing <- setdiff(cols, colnames(meta))
  if (length(missing)) {
    stop("Missing strict-stratum columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invalid <- vapply(meta[, cols, drop = FALSE], function(value) {
    anyNA(value) || any(!nzchar(trimws(as.character(value))))
  }, logical(1))
  if (any(invalid)) {
    stop("Strict-stratum columns contain missing values.", call. = FALSE)
  }
  meta$.rc_stratum_id <- rc_make_stratum_id(meta, cols)
  meta
}
