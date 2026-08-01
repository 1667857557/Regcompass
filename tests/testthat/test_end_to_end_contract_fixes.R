test_that("solver status text is normalized", {
  expect_equal(.rc_lp_status("Optimal"), "optimal")
  expect_equal(.rc_lp_status("infeasible"), "infeasible")
})

test_that("infeasible COMPASS scores are NA rather than zero", {
  penalty <- matrix(
    c(1, 2),
    nrow = 1,
    dimnames = list("R1", c("u1", "u2"))
  )
  feasible <- matrix(
    c(TRUE, FALSE),
    nrow = 1,
    dimnames = dimnames(penalty)
  )
  score <- rc_compass_score_from_penalty(
    penalty,
    feasible
  )
  expect_true(is.na(score["R1", "u2"]))
})
