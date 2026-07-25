.rc_bind_frames_fill <- function(x) {
  x <- x[vapply(
    x,
    function(z) is.data.frame(z) && nrow(z) > 0L,
    logical(1)
  )]
  if (!length(x)) return(data.frame())
  columns <- unique(unlist(lapply(x, colnames), use.names = FALSE))
  x <- lapply(x, function(z) {
    missing <- setdiff(columns, colnames(z))
    for (column in missing) z[[column]] <- NA
    z[, columns, drop = FALSE]
  })
  out <- do.call(rbind, x)
  rownames(out) <- NULL
  out
}

.rc_metacell_logcpm <- function(
    counts, scale_factor = 1e6, library_size = NULL) {
  counts <- .rc_as_dgCMatrix(counts)
  if (!is.numeric(scale_factor) || length(scale_factor) != 1L ||
      !is.finite(scale_factor) || scale_factor <= 0) {
    stop("`scale_factor` must be one positive finite number.", call. = FALSE)
  }
  normalization_scope <- "input_matrix_library_size"
  if (is.null(library_size)) {
    library_size <- Matrix::colSums(counts)
  } else {
    normalization_scope <- "full_transcriptome_library_size_before_gpr_filter"
    if (!is.null(names(library_size))) {
      missing <- setdiff(colnames(counts), names(library_size))
      if (length(missing)) {
        stop(
          "`library_size` is missing metacells: ",
          paste(utils::head(missing, 10L), collapse = ", "),
          call. = FALSE
        )
      }
      library_size <- library_size[colnames(counts)]
    }
  }
  library_size <- as.numeric(library_size)
  if (length(library_size) != ncol(counts) ||
      any(!is.finite(library_size)) || any(library_size <= 0)) {
    stop(
      "`library_size` must contain one positive finite value per metacell.",
      call. = FALSE
    )
  }
  scaled <- counts %*% Matrix::Diagonal(x = scale_factor / library_size)
  answer <- log1p(scaled)
  dimnames(answer) <- dimnames(counts)
  attr(answer, "normalization_scope") <- normalization_scope
  attr(answer, "library_size") <- stats::setNames(
    library_size,
    colnames(counts)
  )
  answer
}
