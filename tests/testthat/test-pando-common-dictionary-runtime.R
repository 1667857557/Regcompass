test_that("condition runtime uses the common-dictionary Pando API", {
  fit_text <- paste(
    deparse(body(RegCompassR:::.rc_fit_condition_grns_by_cell_type)),
    collapse = "\n"
  )
  expect_match(fit_text, "Pando::infer_condition_grn", fixed = TRUE)
  expect_match(fit_text, "adjust_method = \"BH\"", fixed = TRUE)
  expect_match(fit_text, "padj_threshold = 0.05", fixed = TRUE)
  expect_match(fit_text, "rank_action", fixed = TRUE)
  expect_match(fit_text, "min_residual_df", fixed = TRUE)
  expect_false(grepl("engine_control|outer_nfolds|inner_nfolds", fit_text))
})

test_that("execution summaries report the fixed-dictionary engine", {
  diagnostics <- data.frame(
    target = c("G1", "G2", "G3"),
    fit_status = c("ok", "rank_deficient", "insufficient_df"),
    stringsAsFactors = FALSE
  )
  summary <- RegCompassR:::.rc_pando_execution_summary(diagnostics)
  expect_identical(
    summary$fit_engine,
    "two_stage_exact_edge_union_fixed_dictionary_glm"
  )
  expect_false(summary$native_runtime_used)
  expect_false(summary$nested_cv_used)
  expect_identical(summary$targets_total, 3L)
  expect_identical(summary$targets_failed, 1L)
})

test_that("Stage 1 rejects retired condition-engine controls", {
  body_text <- paste(deparse(body(rc_regcompass_step_grn)), collapse = "\n")
  for (argument in c(
    "condition_mix", "outer_nfolds", "inner_nfolds",
    "lambda_selection", "engine_control"
  )) {
    expect_match(body_text, argument, fixed = TRUE)
  }
  expect_match(body_text, "Retired condition-GRN parameter", fixed = TRUE)
})