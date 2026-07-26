# Ensure the candidate table and RNA/ATAC-derived predictor columns have one
# deterministic order before the leakage-resistant fitter constructs matrices.

.rc_fit_multitask_target_cv <- .rc_fit_multitask_target

.rc_fit_multitask_target <- function(
    edges, target, rna, atac, meta, condition_col, args) {
  if (is.data.frame(edges) && nrow(edges)) {
    if (!"edge_id" %in% colnames(edges) || anyDuplicated(edges$edge_id)) {
      stop("Target candidate edges require unique `edge_id` values.",
           call. = FALSE)
    }
    edges <- edges[order(as.character(edges$edge_id)), , drop = FALSE]
    rownames(edges) <- NULL
  }
  .rc_fit_multitask_target_cv(
    edges = edges,
    target = target,
    rna = rna,
    atac = atac,
    meta = meta,
    condition_col = condition_col,
    args = args
  )
}
