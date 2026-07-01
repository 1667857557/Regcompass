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
