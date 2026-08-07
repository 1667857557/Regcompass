.rc_drop_zero_count_atac_features <- function(
    object, atac_assay = "ATAC", context = "ATAC normalization") {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  counts <- .rc_get_assay_counts(object, atac_assay)
  totals <- Matrix::rowSums(counts)
  keep <- is.finite(totals) & totals > 0
  if (!any(keep)) {
    stop(context, " has no peaks with a positive total count.", call. = FALSE)
  }
  diagnostics <- list(
    n_input_peaks = nrow(counts),
    n_zero_count_peaks_excluded = sum(!keep),
    n_retained_peaks = sum(keep),
    zero_count_peak_policy = "exclude_before_tfidf_and_pando"
  )
  if (any(!keep)) {
    object[[atac_assay]] <- subset(
      object[[atac_assay]], features = rownames(counts)[keep]
    )
  }
  list(object = object, diagnostics = diagnostics)
}

.rc_same_matrix_layout <- function(x, y) {
  identical(dim(x), dim(y)) &&
    identical(as.character(rownames(x)), as.character(rownames(y))) &&
    identical(as.character(colnames(x)), as.character(colnames(y)))
}

.rc_align_normalized_assay <- function(object, assay, context) {
  value <- .rc_pando_assay_data(object, assay)
  expected <- as.character(colnames(object))
  observed <- as.character(colnames(value))
  if (anyNA(expected) || any(!nzchar(expected)) || anyDuplicated(expected) ||
      anyNA(observed) || any(!nzchar(observed)) || anyDuplicated(observed)) {
    stop(context, " cell identifiers must be unique and non-empty.",
         call. = FALSE)
  }
  missing <- setdiff(expected, observed)
  extra <- setdiff(observed, expected)
  if (length(missing) || length(extra)) {
    stop(
      context, " normalized assay contains different cells from the Seurat object. ",
      "Missing: ", paste(utils::head(missing, 10L), collapse = ", "),
      "; extra: ", paste(utils::head(extra, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  if (!identical(observed, expected)) {
    object <- .rc_set_assay_matrix(
      object = object,
      assay = assay,
      layer = "data",
      new_data = value[, expected, drop = FALSE]
    )
  }
  object
}

.rc_apply_celltype_shared_tfidf <- function(
    object, celltype_col, atac_assay = "ATAC",
    method = 1, scale.factor = 1e4) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  .rc_validate_celltype_metadata(object@meta.data, celltype_col)
  object <- .rc_prepare_seurat_assays(
    object,
    assays = atac_assay,
    required_layers = "counts",
    optional_layers = "data"
  )
  filtered <- .rc_drop_zero_count_atac_features(
    object, atac_assay, "Cell-type-shared TF-IDF"
  )
  object <- filtered$object
  counts <- .rc_get_assay_counts(object, atac_assay)
  units <- colnames(counts)
  meta <- object@meta.data[
    match(units, rownames(object@meta.data)), , drop = FALSE
  ]
  cell_type <- trimws(as.character(meta[[celltype_col]]))
  groups <- split(units, cell_type)
  local_zero <- vapply(groups, function(group_units) {
    sum(Matrix::rowSums(counts[, group_units, drop = FALSE]) <= 0)
  }, integer(1))

  row_index <- vector("list", length(groups))
  col_index <- vector("list", length(groups))
  values <- vector("list", length(groups))
  names(row_index) <- names(col_index) <- names(values) <- names(groups)
  for (group_name in names(groups)) {
    group_units <- groups[[group_name]]
    group_counts <- counts[, group_units, drop = FALSE]
    keep <- Matrix::rowSums(group_counts) > 0
    if (!any(keep)) next
    normalized <- Signac::RunTFIDF(
      group_counts[keep, , drop = FALSE],
      method = method,
      scale.factor = scale.factor,
      verbose = FALSE
    )
    triplet <- methods::as(.rc_as_sparse(normalized), "TsparseMatrix")
    if (!length(triplet@x)) next
    row_index[[group_name]] <- which(keep)[triplet@i + 1L]
    group_columns <- match(group_units, units)
    col_index[[group_name]] <- group_columns[triplet@j + 1L]
    values[[group_name]] <- as.numeric(triplet@x)
  }

  i <- unlist(row_index, use.names = FALSE)
  j <- unlist(col_index, use.names = FALSE)
  x <- unlist(values, use.names = FALSE)
  if (length(x)) {
    tfidf <- Matrix::sparseMatrix(
      i = i, j = j, x = x,
      dims = dim(counts),
      dimnames = dimnames(counts),
      index1 = TRUE,
      giveCsparse = TRUE
    )
  } else {
    tfidf <- Matrix::Matrix(
      0,
      nrow = nrow(counts),
      ncol = ncol(counts),
      sparse = TRUE,
      dimnames = dimnames(counts)
    )
  }
  tfidf <- .rc_as_sparse(tfidf)
  if (!.rc_same_matrix_layout(tfidf, counts)) {
    stop("Triplet TF-IDF reconstruction changed the ATAC layout.",
         call. = FALSE)
  }
  object <- .rc_set_assay_matrix(
    object = object,
    assay = atac_assay,
    layer = "data",
    new_data = tfidf
  )
  object <- .rc_align_normalized_assay(object, atac_assay, "ATAC")
  object@misc$regcompass_atac_normalization <- c(list(
    method = "Signac_TFIDF",
    assembly = "single_triplet_reconstruction",
    scope = "cell_type_all_available_cells",
    celltype_col = celltype_col,
    idf_reference = "all cells of the same cell type",
    n_units_by_celltype = vapply(groups, length, integer(1)),
    n_zero_count_peaks_by_celltype = local_zero,
    celltype_local_zero_peak_policy =
      "retain_as_zero_without_passing_to_RunTFIDF",
    tfidf_method = method,
    scale_factor = scale.factor
  ), filtered$diagnostics)
  object
}

.rc_validate_condition_celltype_metadata <- function(
    metadata, condition_col = "condition", celltype_col = "cell_type",
    require_multiple_conditions = FALSE) {
  .rc_validate_celltype_metadata(metadata, celltype_col)
  if (is.null(condition_col) || !is.character(condition_col) ||
      length(condition_col) != 1L || is.na(condition_col) ||
      !nzchar(trimws(condition_col)) || !condition_col %in% colnames(metadata)) {
    if (isTRUE(require_multiple_conditions)) {
      stop("A valid `condition_col` is required for condition GRN mode.",
           call. = FALSE)
    }
    return(invisible(TRUE))
  }
  if (identical(condition_col, celltype_col)) {
    stop("`condition_col` and `celltype_col` must name different columns.",
         call. = FALSE)
  }
  value <- as.character(metadata[[condition_col]])
  if (anyNA(value) || any(!nzchar(trimws(value))) ||
      any(value != trimws(value))) {
    stop(
      "Condition labels must be complete, non-empty, and free of surrounding whitespace.",
      call. = FALSE
    )
  }
  if (isTRUE(require_multiple_conditions) && length(unique(value)) < 2L) {
    stop("Condition GRN mode requires at least two condition levels.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.rc_normalize_single_cell_grn_object <- function(
    object, condition_col = NULL, celltype_col = "cell_type",
    rna_assay = "RNA", atac_assay = "ATAC") {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  .rc_validate_condition_celltype_metadata(
    object@meta.data,
    condition_col = condition_col,
    celltype_col = celltype_col,
    require_multiple_conditions = FALSE
  )
  object <- .rc_prepare_seurat_assays(
    object,
    assays = c(rna_assay, atac_assay),
    required_layers = "counts",
    optional_layers = "data"
  )
  object <- Seurat::NormalizeData(object, assay = rna_assay, verbose = FALSE)
  object <- .rc_align_normalized_assay(object, rna_assay, "RNA")
  object <- .rc_apply_celltype_shared_tfidf(
    object, celltype_col = celltype_col, atac_assay = atac_assay
  )
  object@misc$regcompass_grn_normalization <- list(
    rna = "global_single_cell_NormalizeData",
    atac = "cell_type_shared_TFIDF_all_available_cells",
    condition_col = condition_col,
    celltype_col = celltype_col,
    seurat_compatibility = object@misc$regcompass_seurat_compatibility
  )
  object
}

.rc_normalize_condition_metacell_object <- function(
    pooled, rna_assay = "RNA", atac_assay = "ATAC") {
  if (!is.list(pooled) || !inherits(pooled$metacell_object, "Seurat") ||
      !is.data.frame(pooled$metacell_meta)) {
    stop("Stage 2 pooled output lacks its canonical metacell object or metadata.",
         call. = FALSE)
  }
  object <- pooled$metacell_object
  if (!"metacell_id" %in% colnames(pooled$metacell_meta)) {
    stop("Stage 2 metacell metadata lack `metacell_id`.", call. = FALSE)
  }
  expected <- trimws(as.character(pooled$metacell_meta$metacell_id))
  if (anyNA(expected) || any(!nzchar(expected)) || anyDuplicated(expected)) {
    stop("Stage 2 metacell IDs must be unique and non-empty.", call. = FALSE)
  }
  observed <- as.character(colnames(object))
  if (anyNA(observed) || any(!nzchar(observed)) || anyDuplicated(observed) ||
      !setequal(observed, expected)) {
    stop("Stage 2 canonical metacell object and metadata contain different IDs.",
         call. = FALSE)
  }
  if (!identical(observed, expected)) {
    object <- subset(object, cells = expected)
  }
  index <- match(colnames(object), expected)
  object@meta.data <- pooled$metacell_meta[index, , drop = FALSE]
  rownames(object@meta.data) <- colnames(object)

  object <- .rc_prepare_seurat_assays(
    object,
    assays = c(rna_assay, atac_assay),
    required_layers = "counts",
    optional_layers = "data"
  )
  object <- Seurat::NormalizeData(object, assay = rna_assay, verbose = FALSE)
  object <- .rc_align_normalized_assay(object, rna_assay, "RNA")
  .rc_apply_celltype_shared_tfidf(
    object, celltype_col = pooled$celltype_col, atac_assay = atac_assay
  )
}
