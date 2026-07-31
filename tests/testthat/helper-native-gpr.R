.old_reaction_capacity <- function(
    gpr_list, gene_score, promiscuity_mode, and_method, or_method
) {
  gene_score <- as.matrix(gene_score)
  rownames(gene_score) <- tolower(trimws(rownames(gene_score)))
  weights <- rc_promiscuity_weight(gpr_list, mode = promiscuity_mode)
  common <- intersect(rownames(gene_score), names(weights))
  gene_score[common, ] <- sweep(
    gene_score[common, , drop = FALSE], 1L, weights[common], "*"
  )
  value <- lapply(gpr_list, function(rule) {
    vapply(seq_len(ncol(gene_score)), function(pool) {
      score <- gene_score[, pool]
      names(score) <- rownames(gene_score)
      rc_reaction_capacity_one(rule, score, and_method, or_method)
    }, numeric(1))
  })
  out <- do.call(rbind, value)
  dimnames(out) <- list(names(gpr_list), colnames(gene_score))
  out
}
