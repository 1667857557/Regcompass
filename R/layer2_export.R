#' Export RegCompassR Layer 2 results as a long table
#' @export
rc_export_layer2_long <- function(layer2, file = NULL) {
  P <- layer2$L2_compass_like_penalty; S <- layer2$L2_compass_like_score; V <- layer2$L2_vmax_internal; F <- layer2$L2_feasible_flag; ST <- layer2$L2_solver_status
  out <- data.frame(reaction_id = rep(rownames(S), times = ncol(S)), unit_id = rep(colnames(S), each = nrow(S)), L2_compass_like_score = as.vector(S), L2_compass_like_penalty = as.vector(P), L2_vmax_internal = as.vector(V), L2_feasible_flag = as.vector(F), L2_solver_status = as.vector(ST), stringsAsFactors = FALSE)
  if (!is.null(file)) utils::write.table(out, file, sep = "\t", quote = FALSE, row.names = FALSE)
  out
}
