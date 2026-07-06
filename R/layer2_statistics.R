#' Fit per-reaction Layer 2 linear models
#' @export
rc_layer2_lm_by_reaction <- function(score_mat, unit_meta, formula) {
  score_mat <- as.matrix(score_mat)
  fits <- lapply(rownames(score_mat), function(r) {
    dat <- data.frame(L2_score = as.numeric(score_mat[r, unit_meta$unit_id]), unit_meta, stringsAsFactors = FALSE)
    stats::lm(formula, data = dat)
  })
  stats::setNames(fits, rownames(score_mat))
}
