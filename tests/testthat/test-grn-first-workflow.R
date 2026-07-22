test_that("GRN-first defaults are encoded in canonical functions", {
  expect_true(exists("rc_regcompass_step_grn", mode = "function"))
  grn_formals <- paste(deparse(formals(.rc_run_condition_single_cell_grns)$pando_infer_args), collapse = " ")
  expect_match(grn_formals, "peak_cor = 0.01", fixed = TRUE)
  expect_null(eval(formals(rc_regcompass_step_metacells)$sample_col))
})

test_that("metacells are condition-only with gamma 75", {
  body_text <- paste(deparse(body(.rc_make_condition_pooled_metacells)), collapse = "\n")
  expect_match(body_text, "gamma <- 75L", fixed = TRUE)
  expect_match(body_text, 'pooling_scope <- "condition_only"', fixed = TRUE)
  expect_match(body_text, "metacell_grouping = condition_col", fixed = TRUE)
  expect_match(body_text, "Sample balancing is not part", fixed = TRUE)
})

test_that("stepwise meta-modules consume GRN and metacells", {
  f <- formals(rc_regcompass_step_meta_modules)
  expect_true(all(c("grn", "metacells", "gem", "outdir") %in% names(f)))
  expect_false("pando_args" %in% names(f))
})
