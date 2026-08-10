.rc_pando_rna_objects <- function(grn_result) {
  condition_objects <- grn_result$pando_grn_data_by_cell_type %||% list()
  if (!length(condition_objects) &&
      inherits(grn_result$pando_grn_data, "GRNData")) {
    condition_objects <- list(condition_legacy = grn_result$pando_grn_data)
  }
  standard_objects <- grn_result$standard_pando_objects %||% list()
  objects <- c(condition_objects, standard_objects)
  objects <- objects[!vapply(objects, is.null, logical(1))]
  if (!length(objects)) {
    stop(
      "Quantitative RNA aggregation requires the cell-level Pando objects retained by Stage 1.",
      call. = FALSE
    )
  }
  if (anyDuplicated(names(objects))) {
    stop("Pando RNA source object names must be unique.", call. = FALSE)
  }
  invalid <- !vapply(objects, inherits, logical(1), what = "GRNData")
  if (any(invalid)) {
    stop(
      "Quantitative RNA aggregation requires GRNData objects for every routed cell type.",
      call. = FALSE
    )
  }
  objects
}

.rc_single_cell_linear_cpm <- function(
    counts, genes, scale_factor = 1e6) {
  if (is.null(dim(counts)) || is.null(rownames(counts)) ||
      is.null(colnames(counts)) || anyDuplicated(rownames(counts)) ||
      anyDuplicated(colnames(counts))) {
    stop("Single-cell RNA counts require unique gene and cell IDs.",
         call. = FALSE)
  }
  if (!is.numeric(scale_factor) || length(scale_factor) != 1L ||
      !is.finite(scale_factor) || scale_factor <= 0) {
    stop("`scale_factor` must be one positive finite number.", call. = FALSE)
  }
  genes <- unique(tolower(trimws(as.character(genes))))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  feature_key <- tolower(rownames(counts))
  if (anyDuplicated(feature_key)) {
    stop("Duplicated RNA genes after case normalization.", call. = FALSE)
  }
  index <- match(genes, feature_key)
  if (anyNA(index)) {
    missing <- genes[is.na(index)]
    stop(
      "Cell-level Pando RNA counts are missing GPR gene(s); first missing gene: ",
      missing[[1L]], ".",
      call. = FALSE
    )
  }
  library_size <- as.numeric(Matrix::colSums(counts))
  if (any(!is.finite(library_size)) || any(library_size <= 0)) {
    stop("Every cell must have a positive finite RNA library size.",
         call. = FALSE)
  }
  selected <- counts[index, , drop = FALSE]
  rownames(selected) <- genes
  expression <- selected %*% Matrix::Diagonal(
    x = scale_factor / library_size
  )
  dimnames(expression) <- list(genes, colnames(counts))
  list(
    expression = expression,
    library_size = stats::setNames(library_size, colnames(counts)),
    scale_factor = scale_factor
  )
}

.rc_equal_mean_supercell_expression <- function(
    cell_expression, membership, units) {
  if (!is.data.frame(membership) ||
      !all(c("cell_id", "metacell_id") %in% colnames(membership)) ||
      anyNA(membership$cell_id) || anyNA(membership$metacell_id) ||
      any(!nzchar(as.character(membership$cell_id))) ||
      any(!nzchar(as.character(membership$metacell_id))) ||
      anyDuplicated(membership$cell_id)) {
    stop("SuperCell membership must map every cell exactly once.",
         call. = FALSE)
  }
  if (is.null(dim(cell_expression)) || is.null(rownames(cell_expression)) ||
      is.null(colnames(cell_expression)) || anyDuplicated(rownames(cell_expression)) ||
      anyDuplicated(colnames(cell_expression))) {
    stop("Cell-level quantitative RNA requires unique gene and cell IDs.",
         call. = FALSE)
  }
  units <- as.character(units)
  if (anyNA(units) || any(!nzchar(units)) || anyDuplicated(units)) {
    stop("Metacell unit IDs must be unique and non-empty.", call. = FALSE)
  }
  cell_ids <- as.character(membership$cell_id)
  if (!setequal(colnames(cell_expression), cell_ids)) {
    stop(
      "Cell-level quantitative RNA and SuperCell membership must cover exactly the same cells.",
      call. = FALSE
    )
  }
  membership <- membership[
    match(colnames(cell_expression), cell_ids), , drop = FALSE
  ]
  unit_index <- match(as.character(membership$metacell_id), units)
  if (anyNA(unit_index)) {
    stop("SuperCell membership contains metacell IDs absent from Layer 1 units.",
         call. = FALSE)
  }
  n_cells <- tabulate(unit_index, nbins = length(units))
  if (any(n_cells < 1L)) {
    stop("Every Layer 1 metacell must contain at least one member cell.",
         call. = FALSE)
  }
  averaging <- Matrix::sparseMatrix(
    i = seq_len(nrow(membership)),
    j = unit_index,
    x = 1 / n_cells[unit_index],
    dims = c(nrow(membership), length(units)),
    dimnames = list(colnames(cell_expression), units)
  )
  expression <- as.matrix(cell_expression %*% averaging)
  dimnames(expression) <- list(rownames(cell_expression), units)
  list(
    expression = expression,
    n_cells = stats::setNames(as.integer(n_cells), units),
    aggregation = "equal_mean_over_single_cell_linear_cpm",
    library_size_weighted = FALSE
  )
}

.rc_quantitative_supercell_rna <- function(
    grn_result, membership, units, genes, rna_assay,
    scale_factor = 1e6) {
  objects <- .rc_pando_rna_objects(grn_result)
  membership_cells <- as.character(membership$cell_id)
  expression_parts <- list()
  library_parts <- list()
  source_rows <- list()
  for (name in names(objects)) {
    object <- objects[[name]]
    cells <- intersect(membership_cells, colnames(object@data))
    if (!length(cells)) next
    counts <- .rc_get_assay_counts(object@data, rna_assay)
    if (!all(cells %in% colnames(counts))) {
      stop("Pando RNA count columns do not match the stored cell IDs.",
           call. = FALSE)
    }
    normalized <- .rc_single_cell_linear_cpm(
      counts[, cells, drop = FALSE], genes, scale_factor = scale_factor
    )
    expression_parts[[name]] <- normalized$expression
    library_parts[[name]] <- normalized$library_size
    source_rows[[name]] <- data.frame(
      source = name,
      n_cells = length(cells),
      stringsAsFactors = FALSE
    )
  }
  if (!length(expression_parts)) {
    stop("No Stage 1 Pando cells overlap the SuperCell membership.",
         call. = FALSE)
  }
  observed_cells <- unlist(lapply(expression_parts, colnames), use.names = FALSE)
  if (anyDuplicated(observed_cells)) {
    duplicated_cells <- unique(observed_cells[duplicated(observed_cells)])
    stop(
      "A cell occurs in more than one routed Pando RNA source; first duplicated cell: ",
      duplicated_cells[[1L]], ".",
      call. = FALSE
    )
  }
  missing_cells <- setdiff(membership_cells, observed_cells)
  extra_cells <- setdiff(observed_cells, membership_cells)
  if (length(missing_cells) || length(extra_cells)) {
    stop(
      "Stage 1 Pando RNA sources and SuperCell membership are not an exact cell partition; ",
      "missing=", length(missing_cells), ", extra=", length(extra_cells), ".",
      call. = FALSE
    )
  }
  cell_expression <- do.call(cbind, expression_parts)
  cell_expression <- cell_expression[, membership_cells, drop = FALSE]
  averaged <- .rc_equal_mean_supercell_expression(
    cell_expression, membership, units
  )
  library_size <- unlist(library_parts, use.names = TRUE)
  library_size <- library_size[membership_cells]
  if (anyNA(library_size)) {
    stop("Cell-level RNA library-size provenance is incomplete.", call. = FALSE)
  }
  list(
    expression = averaged$expression,
    n_cells = averaged$n_cells,
    cell_library_size = library_size,
    scale_factor = scale_factor,
    aggregation = averaged$aggregation,
    library_size_weighted = averaged$library_size_weighted,
    source_summary = do.call(rbind, source_rows),
    model = "supercell_equal_mean_single_cell_linear_cpm_v1"
  )
}
