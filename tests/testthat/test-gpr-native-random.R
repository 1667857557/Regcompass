test_that("compiled GPR capacities match 48 randomized rule sets", {
  genes <- paste0("G", seq_len(12L))
  and_methods <- c("min", "median", "mean")
  or_methods <- c("max", "sum_sqrtK", "prob_or", "sum")
  modes <- c("none", "sqrt", "linear")
  for (seed in seq_len(48L)) {
    set.seed(seed)
    score <- matrix(
      runif(12L * 7L), nrow = 12L,
      dimnames = list(genes, paste0("P", seq_len(7L)))
    )
    score[sample(length(score), 4L)] <- NA_real_
    gpr <- setNames(lapply(seq_len(10L), function(reaction) {
      n_group <- sample(0:4, 1L)
      if (!n_group) return(list())
      lapply(seq_len(n_group), function(group) {
        sample(c(genes, "missing"), sample(1:4, 1L), replace = TRUE)
      })
    }), paste0("R", seq_len(10L)))
    mode <- modes[(seed - 1L) %% length(modes) + 1L]
    and_method <- and_methods[(seed - 1L) %% length(and_methods) + 1L]
    or_method <- or_methods[(seed - 1L) %% length(or_methods) + 1L]
    reference <- .old_reaction_capacity(
      gpr, score, mode, and_method, or_method
    )
    compiled <- rc_reaction_capacity(
      gpr, score, mode, and_method, or_method
    )
    expect_identical(dimnames(compiled), dimnames(reference), info = seed)
    expect_identical(is.na(compiled), is.na(reference), info = seed)
    expect_equal(compiled, reference, tolerance = 1e-15, info = seed)
  }
})
