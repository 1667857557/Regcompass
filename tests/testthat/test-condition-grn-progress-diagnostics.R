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

test_that("condition fit task reports E-star production and separated inference", {
  body_text <- paste(deparse(body(.rc_condition_ridge_fit_task)), collapse = "\n")
  expect_match(body_text, "verbose = show_progress", fixed = TRUE)
  expect_match(
    body_text,
    "phase=pando_Estar_z025_separate_inference_pipeline",
    fixed = TRUE
  )
  expect_match(
    body_text,
    "Condition-GRN E-star/separate-inference fit failed for cell type",
    fixed = TRUE
  )
  expect_match(
    body_text,
    "phase=pando_Estar_z025_separate_inference_complete",
    fixed = TRUE
  )
  expect_match(body_text, "condition_inference_estimable_rows=", fixed = TRUE)
  expect_match(body_text, "edge_inference_estimable=", fixed = TRUE)
  expect_match(body_text, "regcompass_edges=", fixed = TRUE)
  expect_false(grepl("JSE", body_text, fixed = TRUE))
  expect_false(grepl("condition_significant_rows=", body_text, fixed = TRUE))
})
