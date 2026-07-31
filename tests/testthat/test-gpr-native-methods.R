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

test_that("compiled GPR capacities preserve every aggregation method", {
  gpr <- list(
    R1 = list(c("A", "B")),
    R2 = list(c("C"), c("D", "missing")),
    R3 = list(character()),
    R4 = list(c("A", "A"), c("E"))
  )
  score <- matrix(
    c(0.1, 0.5, 0.9, 0.2, 0.7,
      0.8, 0.4, 0.3, 0.6, NA_real_),
    nrow = 5L,
    dimnames = list(c("A", "B", "C", "D", "E"), c("P1", "P2"))
  )
  for (promiscuity in c("none", "sqrt", "linear")) {
    for (and_method in c("min", "median", "mean")) {
      for (or_method in c("max", "sum_sqrtK", "prob_or", "sum")) {
        reference <- .old_reaction_capacity(
          gpr, score, promiscuity, and_method, or_method
        )
        compiled <- rc_reaction_capacity(
          gpr, score, promiscuity, and_method, or_method
        )
        expect_identical(dimnames(compiled), dimnames(reference))
        expect_identical(is.na(compiled), is.na(reference))
        expect_equal(compiled, reference, tolerance = 1e-15)
      }
    }
  }
})
