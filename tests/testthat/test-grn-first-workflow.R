test_that("GRN-first defaults are encoded in canonical functions", {
  expect_true(exists("rc_regcompass_step_grn", mode = "function"))
  grn_formals <- paste(
    deparse(formals(.rc_fit_condition_grns_by_cell_type)$pando_infer_args),
    collapse = " "
  )
  expect_match(grn_formals, "peak_cor = 0", fixed = TRUE)
  expect_null(eval(formals(rc_regcompass_step_metacells)$sample_col))
})

test_that("metacells are condition-by-cell-type with gamma 30", {
  body_text <- paste(deparse(body(.rc_make_condition_celltype_metacells)), collapse = "\n")
  expect_match(body_text, "gamma <- 30L", fixed = TRUE)
  expect_identical(eval(formals(.rc_build_supercell2_strata)$gamma), 30)
  expect_identical(eval(formals(rc_make_supercell2_metacells)$gamma), 30)
  expect_match(
    body_text, 'pooling_scope <- "condition_by_cell_type"', fixed = TRUE
  )
  expect_match(
    body_text, "metacell_grouping = c(condition_col, celltype_col)",
    fixed = TRUE
  )
  expect_match(body_text, "supercell_stratum_col", fixed = TRUE)
})

test_that("stepwise meta-modules consume GRN and metacells", {
  f <- formals(rc_regcompass_step_meta_modules)
  expect_true(all(c("grn", "metacells", "gem", "outdir") %in% names(f)))
  expect_false("pando_args" %in% names(f))
})
