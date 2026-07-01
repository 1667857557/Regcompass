#' Minimal pool diagnostics for the simplified Layer 1 workflow
#' @export
rc_pool_diagnostics <- function(pool_map,
                                rna_counts = NULL,
                                sample_col = "sample_id",
                                celltype_col = "cell_type",
                                condition_col = NULL,
                                gpr_genes = NULL,
                                BPPARAM = NULL) {
  if (!is.data.frame(pool_map)) stop("`pool_map` must be a data.frame.", call. = FALSE)
  missing_cols <- setdiff(c("pool_id", "cell_id"), colnames(pool_map))
  if (length(missing_cols) > 0L) stop("`pool_map` is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)

  active <- pool_map
  if ("skipped" %in% colnames(active)) active <- active[!active$skipped, , drop = FALSE]
  active <- active[!is.na(active$pool_id), , drop = FALSE]
  if (!is.null(rna_counts)) rc_validate_diagnostic_matrix(rna_counts, active$cell_id, "rna_counts")

  rna_depth <- if (is.null(rna_counts)) NULL else Matrix::colSums(rna_counts)
  gpr_genes <- rc_match_matrix_features(gpr_genes, rna_counts)
  pool_ids <- unique(active$pool_id)

  pieces <- rc_pool_lapply(pool_ids, function(pid) {
    one <- active[active$pool_id == pid, , drop = FALSE]
    cells <- one$cell_id
    data.frame(
      pool_id = pid,
      sample_id = rc_pool_unique_value(one, sample_col),
      cell_type = rc_pool_unique_value(one, celltype_col),
      condition = rc_pool_unique_value(one, condition_col),
      n_cells = length(cells),
      low_power_pool = if ("low_power_pool" %in% colnames(one)) any(one$low_power_pool) else NA,
      no_within_group_pool_replicate = if ("no_within_group_pool_replicate" %in% colnames(one)) any(one$no_within_group_pool_replicate) else NA,
      RNA_depth_mean = if (is.null(rna_depth)) NA_real_ else mean(rna_depth[cells], na.rm = TRUE),
      GPR_gene_detection_rate = rc_feature_detection_mean(rna_counts, gpr_genes, cells),
      stringsAsFactors = FALSE
    )
  }, BPPARAM = BPPARAM)

  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

rc_validate_diagnostic_matrix <- function(mat, cells, name) {
  if (is.null(dim(mat)) || length(dim(mat)) != 2L) stop("`", name, "` must be a two-dimensional feature-by-cell matrix.", call. = FALSE)
  if (is.null(colnames(mat))) stop("`", name, "` must have cell IDs in colnames().", call. = FALSE)
  missing_cells <- setdiff(cells, colnames(mat))
  if (length(missing_cells) > 0L) stop("Some pool_map cell IDs are absent from `", name, "`: ", paste(utils::head(missing_cells), collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

rc_pool_unique_value <- function(df, col) {
  if (is.null(col) || is.na(col) || !nzchar(col) || !col %in% colnames(df)) return(NA_character_)
  vals <- unique(as.character(df[[col]]))
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0L) return(NA_character_)
  if (length(vals) > 1L) return(paste(vals, collapse = ";"))
  vals[[1]]
}

rc_match_matrix_features <- function(features, mat) {
  if (is.null(features) || is.null(mat) || is.null(rownames(mat))) return(character(0))
  lower_map <- stats::setNames(rownames(mat), tolower(rownames(mat)))
  matched <- lower_map[tolower(unique(features))]
  unique(stats::na.omit(as.character(matched)))
}

rc_feature_detection_mean <- function(mat, features, cells) {
  if (is.null(mat) || length(features) == 0L || length(cells) == 0L) return(NA_real_)
  mean(mat[features, cells, drop = FALSE] > 0, na.rm = TRUE)
}
