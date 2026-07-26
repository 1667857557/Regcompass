test_that("GRN-first defaults are encoded in canonical functions", {
  expect_true(exists("rc_regcompass_step_grn", mode = "function"))
  grn_formals <- paste(
    deparse(formals(.rc_run_condition_single_cell_grns)$pando_infer_args),
    collapse = " "
  )
  expect_match(grn_formals, "peak_cor = 0.01", fixed = TRUE)
  expect_false("sample_col" %in% names(formals(rc_regcompass_step_metacells)))
})

test_that("current SuperCell2 builder uses strata and label without sample adapters", {
  builder <- formals(rc_make_supercell2_metacells)
  expect_true(all(c("strata_cols", "label_col") %in% names(builder)))
  expect_false(any(c("sample_col", "pool_col") %in% names(builder)))
  expect_identical(eval(builder$strata_cols), "condition")
  expect_identical(builder$gamma, 30)

  body_text <- paste(
    deparse(body(.rc_make_condition_pooled_metacells)), collapse = "\n"
  )
  expect_match(body_text, "strata_cols = condition_col", fixed = TRUE)
  expect_match(body_text, "label_col = celltype_col", fixed = TRUE)
  expect_match(body_text, 'pooling_scope <- "condition_only"', fixed = TRUE)
  expect_false(grepl("condition_pool_id", body_text, fixed = TRUE))
  expect_false(grepl("sample_col", body_text, fixed = TRUE))
})

test_that("stepwise meta-modules consume GRN and metacells", {
  f <- formals(rc_regcompass_step_meta_modules)
  expect_true(all(c("grn", "metacells", "gem", "outdir") %in% names(f)))
  expect_false("pando_args" %in% names(f))
})
