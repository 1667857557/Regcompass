#' Export a reaction-by-sample matrix
#' @export
rc_export_sample_matrix <- function(sample_matrix, file) {
  sample_matrix <- as.matrix(sample_matrix)
  out <- data.frame(reaction_id = rownames(sample_matrix), sample_matrix, check.names = FALSE)
  utils::write.table(out, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(file)
}

#' Export a pool- or sample-level long table from a reaction-by-column matrix
#' @export
rc_export_long_table <- function(score_mat, file, value_col = "value") {
  score_mat <- as.matrix(score_mat)
  long <- data.frame(
    reaction_id = rep(rownames(score_mat), times = ncol(score_mat)),
    column_id = rep(colnames(score_mat), each = nrow(score_mat)),
    value = as.vector(score_mat),
    stringsAsFactors = FALSE
  )
  names(long)[3] <- value_col
  utils::write.table(long, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(long)
}

#' Sample-level long summary with pool/cell-count diagnostics
rc_sample_summary <- function(score_mat, pool_meta, sample_col = "sample_id", celltype_col = "cell_type", condition_col = NULL) {
  score_mat <- as.matrix(score_mat)
  required <- c("pool_id", sample_col, celltype_col)
  missing <- setdiff(required, colnames(pool_meta))
  if (length(missing) > 0L) stop("`pool_meta` is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  pool_meta <- pool_meta[pool_meta$pool_id %in% colnames(score_mat), , drop = FALSE]
  group_cols <- c(sample_col, condition_col, celltype_col)
  group_cols <- group_cols[!is.null(group_cols) & !is.na(group_cols)]
  keys <- interaction(pool_meta[, group_cols, drop = FALSE], drop = TRUE, sep = "|")
  pieces <- lapply(split(pool_meta, keys), function(pm) {
    vals <- score_mat[, pm$pool_id, drop = FALSE]
    data.frame(
      group_id = unique(keys[match(pm$pool_id, pool_meta$pool_id)])[1],
      reaction_id = rownames(score_mat),
      value = matrixStats::rowMedians(vals, na.rm = TRUE),
      n_pools_used = ncol(vals),
      n_cells_used = if ("n_cells" %in% colnames(pm)) sum(pm$n_cells, na.rm = TRUE) else NA_real_,
      low_power_group_flag = if ("low_power_pool" %in% colnames(pm)) any(pm$low_power_pool, na.rm = TRUE) else NA,
      single_pool_group_flag = ncol(vals) == 1L,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, pieces)
}
