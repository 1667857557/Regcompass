test_that("condition target reliability uses sqrt R-squared after final edge gate", {
  fit <- list(
    fit = data.frame(
      target = c("G1", "G1", "G2", "G2", "G3", "G4"),
      condition = c("A", "B", "A", "B", "A", "A"),
      rsq = c(0.25, 0.81, 0.04, NA, 1.44, -0.2),
      fit_status = c("ok", "ok", "rank_deficient", "ok", "ok", "ok"),
      stringsAsFactors = FALSE
    ),
    coefficients = data.frame(
      target = c("G1", "G1", "G2", "G2", "G3", "G4"),
      condition = c("A", "B", "A", "B", "A", "A"),
      significant = c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
      stringsAsFactors = FALSE
    )
  )

  reliability <- .rc_condition_target_reliability(fit)

  expect_equal(reliability$reliability, c(0.5, 0.9, NA, NA, 1, 0))
  expect_equal(
    reliability$n_significant_edges,
    c(1L, 1L, 1L, 0L, 1L, 1L)
  )
})

test_that("condition target reliability counts multiple active edges once at target level", {
  fit <- list(
    fit = data.frame(
      target = "G1", condition = "A", rsq = 0.49,
      fit_status = "ok", stringsAsFactors = FALSE
    ),
    coefficients = data.frame(
      target = c("G1", "g1", "G1"),
      condition = c("A", "A", "A"),
      significant = c(TRUE, TRUE, FALSE),
      stringsAsFactors = FALSE
    )
  )

  reliability <- .rc_condition_target_reliability(fit)

  expect_equal(reliability$n_significant_edges, 2L)
  expect_equal(reliability$reliability, 0.7)
})
