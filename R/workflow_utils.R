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

.rc_validate_celltype_metadata <- function(
    metadata, celltype_col = "cell_type") {
  if (!is.data.frame(metadata)) {
    stop("Cell metadata must be a data frame.", call. = FALSE)
  }
  if (!is.character(celltype_col) || length(celltype_col) != 1L ||
      is.na(celltype_col) || !nzchar(trimws(celltype_col)) ||
      !celltype_col %in% colnames(metadata)) {
    stop("`celltype_col` must identify one metadata column.", call. = FALSE)
  }
  value <- as.character(metadata[[celltype_col]])
  if (anyNA(value) || any(!nzchar(trimws(value))) ||
      any(value != trimws(value))) {
    stop(
      "Cell-type labels must be complete, non-empty, and free of surrounding whitespace.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.rc_internal_condition_col <- function(metadata) {
  name <- ".regcompass_condition"
  while (name %in% colnames(metadata)) name <- paste0(name, "_")
  name
}

.rc_resolve_condition_design <- function(
    object, condition_col = "condition", constant_label = "all") {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  metadata <- object@meta.data
  omitted <- is.null(condition_col) || !is.character(condition_col) ||
    length(condition_col) != 1L || is.na(condition_col) ||
    !nzchar(trimws(condition_col))
  if (omitted) {
    effective <- .rc_internal_condition_col(metadata)
    object@meta.data[[effective]] <- rep(constant_label, nrow(metadata))
    return(list(
      object = object,
      requested_condition_col = condition_col,
      condition_col = effective,
      condition_levels = constant_label,
      analysis_mode = "standard_pando",
      condition_supplied = FALSE,
      fallback_reason = "condition_col_not_supplied"
    ))
  }
  condition_col <- trimws(condition_col)
  if (!condition_col %in% colnames(metadata)) {
    stop(
      "Explicitly requested `condition_col` was not found: `",
      condition_col, "`.",
      call. = FALSE
    )
  }
  values <- as.character(metadata[[condition_col]])
  if (anyNA(values) || any(!nzchar(trimws(values))) ||
      any(values != trimws(values))) {
    stop(
      "Condition labels must be complete, non-empty, and free of surrounding whitespace.",
      call. = FALSE
    )
  }
  levels <- unique(values)
  list(
    object = object,
    requested_condition_col = condition_col,
    condition_col = condition_col,
    condition_levels = levels,
    analysis_mode = if (length(levels) >= 2L) {
      "condition_grn"
    } else {
      "standard_pando"
    },
    condition_supplied = TRUE,
    fallback_reason = if (length(levels) >= 2L) {
      NA_character_
    } else {
      "fewer_than_two_condition_levels"
    }
  )
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
