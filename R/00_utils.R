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
.rc_as_dgCMatrix <- function(x) {
  if (inherits(x, "dgCMatrix")) return(x)
  if (inherits(x, "sparseMatrix")) {
    out <- methods::as(methods::as(x, "generalMatrix"), "CsparseMatrix")
    if (!inherits(out, "dgCMatrix")) {
      out <- out * 1
      out <- methods::as(methods::as(out, "generalMatrix"), "CsparseMatrix")
    }
    if (!inherits(out, "dgCMatrix")) {
      stop("Could not coerce sparse input to `dgCMatrix`.", call. = FALSE)
    }
    return(out)
  }
  dense <- as.matrix(x)
  storage.mode(dense) <- "double"
  out <- Matrix::Matrix(dense, sparse = TRUE)
  if (!inherits(out, "dgCMatrix")) {
    out <- methods::as(methods::as(out, "generalMatrix"), "CsparseMatrix")
  }
  out
}

.rc_as_sparse <- function(x) {
  .rc_as_dgCMatrix(x)
}

# Clamp finite numerical evidence to the unit interval.
