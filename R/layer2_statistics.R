#' Fit per-reaction Layer 2 linear models
#' @export
rc_layer2_lm_by_reaction <- function(
    score_mat, unit_meta, formula,
    sample_col = "sample_id", celltype_col = "cell_type",
    condition_col = "condition") {
  score_mat <- as.matrix(score_mat)
  if (is.null(rownames(score_mat)) || is.null(colnames(score_mat))) {
    stop("`score_mat` must have reaction rows and unit columns.", call. = FALSE)
  }
  if (!is.data.frame(unit_meta) || !"unit_id" %in% colnames(unit_meta)) {
    stop("`unit_meta` must contain `unit_id`.", call. = FALSE)
  }
  unit_meta <- unit_meta[match(colnames(score_mat), as.character(unit_meta$unit_id)), , drop = FALSE]
  if (anyNA(unit_meta$unit_id)) stop("`unit_meta` is missing score columns.", call. = FALSE)
  required <- c(sample_col, condition_col)
  missing <- setdiff(required, colnames(unit_meta))
  if (length(missing)) stop("`unit_meta` is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  celltypes <- if (celltype_col %in% colnames(unit_meta)) unique(as.character(unit_meta[[celltype_col]])) else "all"
  fits <- list()
  for (cell_type in celltypes) {
    columns <- if (identical(cell_type, "all")) seq_len(ncol(score_mat)) else
      which(as.character(unit_meta[[celltype_col]]) == cell_type)
    aggregated <- .rc_aggregate_microcompass_samples(
      score_mat[, columns, drop = FALSE], unit_meta[columns, , drop = FALSE],
      sample_col = sample_col, condition_col = condition_col
    )
    for (reaction in rownames(aggregated$score)) {
      data <- aggregated$meta
      data$L2_score <- as.numeric(aggregated$score[reaction, rownames(data)])
      fits[[paste(cell_type, reaction, sep = "::")]] <- stats::lm(formula, data = data)
    }
  }
  fits
}
