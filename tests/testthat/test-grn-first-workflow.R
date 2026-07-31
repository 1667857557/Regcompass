test_that("GRN-first defaults are encoded in canonical functions", {
  expect_true(exists("rc_regcompass_step_grn", mode = "function"))
  grn_formals <- paste(
    deparse(formals(.rc_fit_condition_grns_by_cell_type)$pando_infer_args),
    collapse = " "
  )
  expect_match(grn_formals, "peak_cor = 0", fixed = TRUE)
  expect_true("condition_col" %in% names(formals(rc_regcompass_step_metacells)))
  expect_false("sample_col" %in% names(formals(rc_regcompass_step_metacells)))
})

test_that("metacells use grouped WNN and post-clustering condition splitting", {
  defaults <- .rc_condition_metacell_defaults()
  expect_identical(defaults$gamma, 20L)
  build_text <- paste(
    deparse(body(.rc_build_grouped_wnn_membership)), collapse = "\n"
  )
  expect_match(build_text, "SCimplify_by_graph_group", fixed = TRUE)
  expect_match(build_text, "cell.graph.group", fixed = TRUE)
  expect_match(build_text, "cell.split.condition", fixed = TRUE)
  expect_match(build_text, 'assay = c(rna_assay, atac_assay)', fixed = TRUE)
  expect_false(exists(".rc_native_supercell_membership", inherits = TRUE))
})

test_that("stepwise meta-modules consume GRN and metacells", {
  f <- formals(rc_regcompass_step_meta_modules)
  expect_true(all(c("grn", "metacells", "gem", "outdir") %in% names(f)))
  expect_false("pando_args" %in% names(f))
})
