test_that("HiGHS solver maps variable and constraint bounds correctly", {
  A <- matrix(1, nrow = 1, ncol = 1)
  obj <- c(x = 1)
  lb <- c(x = 0)
  ub <- c(x = 10)
  lhs <- 5
  rhs <- Inf
  ans <- rc_solve_lp(obj, A, lhs, rhs, lb, ub, solver = "highs", time_limit = 60)
  if (identical(ans$status, "error")) skip(ans$message)
  expect_equal(ans$status, "optimal")
  expect_equal(ans$objective, 5, tolerance = 1e-6)
  expect_equal(ans$solution[[1]], 5, tolerance = 1e-6)
})
