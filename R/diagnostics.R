#' Compute pool-level diagnostics for sample-aware micropools
#'
#' @param pool_map A pool assignment data.frame from [rc_make_pools()]. Must
#' contain `pool_id` and `cell_id`; grouping columns are reused when present.
#' @param rna_counts Optional raw RNA count matrix with genes/features in rows and
#' cells in columns.
#' @param atac_counts Optional raw ATAC count matrix with peaks/features in rows
#' and cells in columns.
#' @param sample_col Column in `pool_map` containing sample IDs.
#' @param condition_col Optional column in `pool_map` containing condition labels.
#' @param celltype_col Column in `pool_map` containing cell type labels.
#' @param state_col Optional column in `pool_map` containing local state labels.
#' @param metabolic_genes Optional genes used to compute metabolic gene detection.
#' @param gpr_genes Optional genes used to compute GPR gene detection.
#' @param BPPARAM Optional `BiocParallelParam` for pool-wise parallel diagnostics.
#'
#' @return A data.frame with one row per pool and v0.4 diagnostic fields.
#' @export
rc_pool_diagnostics <- function(pool_map,
                                rna_counts = NULL,
                                atac_counts = NULL,
                                sample_col = "sample_id",
                                condition_col = "condition",
                                celltype_col = "cell_type",
                                state_col = NULL,
                                metabolic_genes = NULL,
                                gpr_genes = NULL,
                                BPPARAM = NULL) {
  if (!is.data.frame(pool_map)) {
    stop("`pool_map` must be a data.frame.", call. = FALSE)
  }
  missing_cols <- setdiff(c("pool_id", "cell_id"), colnames(pool_map))
  if (length(missing_cols) > 0) {
    stop("`pool_map` is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  active_pool_map <- pool_map
  if ("skipped" %in% colnames(active_pool_map)) {
    active_pool_map <- active_pool_map[!active_pool_map$skipped, , drop = FALSE]
  }
  active_pool_map <- active_pool_map[!is.na(active_pool_map$pool_id), , drop = FALSE]

  if (!is.null(rna_counts)) {
    rc_validate_diagnostic_matrix(rna_counts, active_pool_map$cell_id, "rna_counts")
  }
  if (!is.null(atac_counts)) {
    rc_validate_diagnostic_matrix(atac_counts, active_pool_map$cell_id, "atac_counts")
  }

  rna_depth <- if (is.null(rna_counts)) NULL else Matrix::colSums(rna_counts)
  atac_depth <- if (is.null(atac_counts)) NULL else Matrix::colSums(atac_counts)
  metabolic_genes <- rc_match_matrix_features(metabolic_genes, rna_counts)
  gpr_genes <- rc_match_matrix_features(gpr_genes, rna_counts)

  pool_ids <- unique(active_pool_map$pool_id)
  pieces <- rc_parallel_lapply(pool_ids, function(pid) {
    rows <- active_pool_map$pool_id == pid
    cells <- active_pool_map$cell_id[rows]
    one <- active_pool_map[rows, , drop = FALSE]

    data.frame(
      pool_id = pid,
      sample_id = rc_pool_unique_value(one, sample_col),
      condition = rc_pool_unique_value(one, condition_col),
      cell_type = rc_pool_unique_value(one, celltype_col),
      local_state = rc_pool_unique_value(one, state_col),
      n_cells = length(cells),
      low_power_pool = if ("low_power_pool" %in% colnames(one)) any(one$low_power_pool) else NA,
      RNA_depth_mean = if (is.null(rna_depth)) NA_real_ else mean(rna_depth[cells], na.rm = TRUE),
      ATAC_depth_mean = if (is.null(atac_depth)) NA_real_ else mean(atac_depth[cells], na.rm = TRUE),
      metabolic_gene_detection_rate = rc_feature_detection_mean(rna_counts, metabolic_genes, cells),
      GPR_gene_detection_rate = rc_feature_detection_mean(rna_counts, gpr_genes, cells),
      stringsAsFactors = FALSE
    )
  }, BPPARAM = BPPARAM)

  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

rc_validate_diagnostic_matrix <- function(mat, cells, name) {
  if (is.null(dim(mat)) || length(dim(mat)) != 2L) {
    stop("`", name, "` must be a two-dimensional feature-by-cell matrix.", call. = FALSE)
  }
  if (is.null(colnames(mat))) {
    stop("`", name, "` must have cell IDs in colnames().", call. = FALSE)
  }
  missing_cells <- setdiff(cells, colnames(mat))
  if (length(missing_cells) > 0) {
    stop("Some pool_map cell IDs are absent from `", name, "`: ", paste(utils::head(missing_cells), collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

rc_pool_unique_value <- function(df, col) {
  if (is.null(col) || is.na(col) || !nzchar(col) || !col %in% colnames(df)) {
    return(NA_character_)
  }
  vals <- unique(as.character(df[[col]]))
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0L) return(NA_character_)
  if (length(vals) > 1L) return(paste(vals, collapse = ";"))
  vals[[1]]
}

rc_match_matrix_features <- function(features, mat) {
  if (is.null(features) || is.null(mat) || is.null(rownames(mat))) {
    return(character(0))
  }
  lower_map <- stats::setNames(rownames(mat), tolower(rownames(mat)))
  matched <- lower_map[tolower(unique(features))]
  unique(stats::na.omit(as.character(matched)))
}

rc_feature_detection_mean <- function(mat, features, cells) {
  if (is.null(mat) || length(features) == 0L || length(cells) == 0L) {
    return(NA_real_)
  }
  mean(mat[features, cells, drop = FALSE] > 0, na.rm = TRUE)
}
