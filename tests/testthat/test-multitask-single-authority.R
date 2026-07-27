test_that("multitask defaults and final fitter have one authoritative contract", {
  defaults_body <- paste(deparse(body(.rc_multitask_grn_defaults)), collapse = "\n")
  fitter_body <- paste(deparse(body(.rc_fit_multitask_target)), collapse = "\n")
  utility_body <- paste(deparse(body(.rc_residualize_matrix)), collapse = "\n")

  expect_match(defaults_body, "n_bootstrap = 100L", fixed = TRUE)
  expect_match(fitter_body, ".rc_fit_multitask_target_direct", fixed = TRUE)
  expect_match(fitter_body, ".rc_cv_predictive_gate", fixed = TRUE)
  expect_false(grepl("global_penalty_factor", utility_body, fixed = TRUE))

  namespace <- asNamespace("RegCompassR")
  for (obsolete in c(
    ".rc_fit_multitask_target_cv",
    ".rc_fit_multitask_target_core",
    ".rc_fit_multitask_target_pre_policy",
    ".rc_multitask_grn_defaults_core",
    ".rc_validate_multitask_grn_args_core"
  )) {
    expect_false(exists(obsolete, envir = namespace, inherits = FALSE))
  }
})

test_that("canonical peak-to-gene policy is GREAT", {
  defaults <- .rc_validate_canonical_pando_design_args()
  expect_identical(defaults$peak_to_gene_method, "GREAT")
  expect_identical(defaults$extend, 1000000)
})
