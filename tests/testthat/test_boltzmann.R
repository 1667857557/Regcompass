test_that("rc_boltzmann_minavg is bounded by min and mean", {
  scores <- c(0.2, 0.5, 0.9)
  out <- rc_boltzmann_minavg(scores, tau = 0.08)
  expect_gte(out, min(scores))
  expect_lte(out, mean(scores))
})

test_that("rc_reaction_capacity_one handles AND and OR semantics", {
  parsed <- list(c("g1", "g2"), "g3")
  gene_score <- c(g1 = 0.2, g2 = 0.8, g3 = 0.5)
  and_part <- rc_boltzmann_minavg(c(0.2, 0.8), tau = 0.08)
  expect_equal(rc_reaction_capacity_one(parsed, gene_score, tau = 0.08), and_part + 0.5)
})
