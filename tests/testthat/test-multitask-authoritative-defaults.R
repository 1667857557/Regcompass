test_that("authoritative multitask defaults match the documented policy", {
  defaults <- .rc_multitask_grn_defaults()
  expect_identical(defaults$deviation_penalty_factor, 1)
  expect_identical(defaults$n_bootstrap, 100L)
  expect_identical(defaults$min_bootstrap_success_fraction, 0.8)
})
