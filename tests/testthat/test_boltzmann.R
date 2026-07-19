test_that("normalized Boltzmann soft-min is bounded by min and mean", {
  scores <- c(0.2, 0.5, 0.9)
  out <- rc_boltzmann_minavg(scores, tau = 0.08)
  expect_gte(out, min(scores))
  expect_lte(out, mean(scores))
})

test_that("normalized Boltzmann soft-min is monotone in every subunit", {
  scores <- c(0.2, 0.5, 0.9)
  baseline <- rc_boltzmann_minavg(scores, tau = 0.20)
  for (index in seq_along(scores)) {
    increased <- scores
    increased[[index]] <- increased[[index]] + 0.1
    expect_gte(
      rc_boltzmann_minavg(increased, tau = 0.20),
      baseline
    )
  }
  expect_equal(
    rc_boltzmann_minavg(rep(0.4, 3), tau = 0.20),
    0.4
  )
})

test_that("rc_reaction_capacity_one handles AND and OR semantics", {
  parsed <- list(c("g1", "g2"), "g3")
  gene_score <- c(g1 = 0.2, g2 = 0.8, g3 = 0.5)
  and_part <- rc_boltzmann_minavg(c(0.2, 0.8), tau = 0.08)
  expect_equal(
    rc_reaction_capacity_one(
      parsed, gene_score, tau = 0.08, and_method = "boltzmann",
      or_method = "sum_sqrtK"
    ),
    (and_part + 0.5) / sqrt(2)
  )
  expect_equal(
    rc_reaction_capacity_one(parsed, gene_score, tau = 0.08),
    max(and_part, 0.5)
  )
})

test_that("AND aggregation supports min, Boltzmann soft-min, and mean", {
  scores <- c(0.2, 0.5, 0.9)
  expect_equal(rc_and_capacity(scores, method = "min"), min(scores))
  expect_equal(rc_and_capacity(scores, method = "mean"), mean(scores))
  boltz <- rc_and_capacity(scores, method = "boltzmann", tau = 0.20)
  expect_gte(boltz, min(scores))
  expect_lte(boltz, mean(scores))
  expect_equal(
    rc_or_capacity(c(0.2, 0.5, NA), method = "sum_sqrtK"),
    0.7 / sqrt(2)
  )
  expect_equal(rc_or_capacity(c(0.2, 0.5, NA)), 0.5)
})
