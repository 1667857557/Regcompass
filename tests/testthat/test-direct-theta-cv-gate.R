test_that("final multitask fitter applies the strict predictive gate", {
  body_text <- paste(deparse(body(.rc_fit_multitask_target)), collapse = "\n")
  expect_match(body_text, ".rc_fit_multitask_target_direct", fixed = TRUE)
  expect_match(body_text, ".rc_cv_predictive_gate", fixed = TRUE)
})
