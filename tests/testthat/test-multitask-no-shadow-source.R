test_that("loaded multitask contract exposes no obsolete shadow aliases", {
  namespace <- asNamespace("RegCompassR")
  obsolete <- c(
    ".rc_fit_multitask_target_cv",
    ".rc_fit_multitask_target_core",
    ".rc_fit_multitask_target_pre_policy",
    ".rc_multitask_grn_defaults_core",
    ".rc_validate_multitask_grn_args_core"
  )
  expect_false(any(vapply(
    obsolete,
    exists,
    logical(1),
    envir = namespace,
    inherits = FALSE
  )))
})
