test_that("condition target pool disables generic BiocParallel progress bars", {
  skip_if_not_installed("BiocParallel")
  param <- BiocParallel::SnowParam(
    workers = 2L,
    type = "SOCK",
    progressbar = TRUE,
    exportglobals = TRUE,
    exportvariables = TRUE
  )
  attr(param, "regcompass_worker_limit") <- 2L
  quiet <- .rc_pando_quiet_target_bpparam(param)
  expect_true(methods::is(quiet, "SnowParam"))
  expect_equal(BiocParallel::bpnworkers(quiet), 2L)
  expect_false(BiocParallel::bpprogressbar(quiet))
  expect_equal(.rc_bpparam_worker_limit(quiet), 2L)
})

test_that("condition fit task exposes semantic Pando progress", {
  body_text <- paste(deparse(body(.rc_condition_ridge_fit_task)), collapse = "\n")
  expect_match(body_text, "verbose = show_progress", fixed = TRUE)
  expect_match(body_text, "phase=pando_condition_pipeline", fixed = TRUE)
  expect_match(body_text, "Condition-GRN ridge fit failed for cell type", fixed = TRUE)
  expect_match(body_text, "pando_condition_pipeline_complete", fixed = TRUE)
  expect_match(body_text, "active_edges=", fixed = TRUE)
})
