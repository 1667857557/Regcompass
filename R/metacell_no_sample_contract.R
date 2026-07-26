# Public pool-based adapter loaded immediately after metacell.R.

.rc_make_supercell2_metacells_raw <- rc_make_supercell2_metacells

.rc_remove_sample_fields <- function(x) {
  if (is.data.frame(x)) {
    keep <- !grepl("sample", colnames(x), ignore.case = TRUE)
    return(x[, keep, drop = FALSE])
  }
  if (!is.list(x)) return(x)
  names_to_remove <- names(x)[grepl("sample", names(x), ignore.case = TRUE)]
  x[names_to_remove] <- NULL
  for (name in names(x)) {
    if (is.data.frame(x[[name]])) {
      x[[name]] <- .rc_remove_sample_fields(x[[name]])
    }
  }
  x
}

#' Build SuperCell2 metacells from explicit pool strata
#'
#' This lower-level builder uses a generic pool identifier rather than a
#' biological-sample column. The canonical RegCompass workflow supplies an
#' internal condition-derived pool and performs no sample balancing or sample
#' composition analysis.
#'
#' @param object Paired RNA+ATAC Seurat object.
#' @param outdir Output directory.
#' @param pool_col Metadata column defining the leading technical pool stratum.
#' @param condition_col Condition metadata column.
#' @param celltype_col Cell-type stratum column used by the low-level builder.
#' @param label_col Optional label supplied to SuperCell2 before aggregation.
#' @param ... Remaining non-sample metacell arguments accepted by the legacy
#'   computational backend.
#' @export
rc_make_supercell2_metacells <- function(
    object, outdir, pool_col,
    condition_col = "condition",
    celltype_col = "cell_type",
    label_col = NULL,
    ...) {
  if (!is.character(pool_col) || length(pool_col) != 1L ||
      is.na(pool_col) || !nzchar(pool_col)) {
    stop("`pool_col` must name one complete metadata column.", call. = FALSE)
  }
  dots <- list(...)
  forbidden <- names(dots)[grepl("sample", names(dots), ignore.case = TRUE)]
  if (length(forbidden)) {
    stop(
      "Sample-specific metacell arguments have been removed: ",
      paste(forbidden, collapse = ", "), call. = FALSE
    )
  }
  call_args <- c(
    list(object, outdir, pool_col),
    list(
      condition_col = condition_col,
      celltype_col = celltype_col,
      label_col = label_col,
      sample_balance = FALSE
    ),
    dots
  )
  answer <- do.call(.rc_make_supercell2_metacells_raw, call_args)
  answer <- .rc_remove_sample_fields(answer)
  answer$pool_col <- pool_col
  answer$pooling_semantics <- "generic_pool_without_sample_balancing"
  answer
}
