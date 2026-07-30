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

test_that("metacells use cell-type graph groups and condition splitting", {
  body_text <- paste(
    deparse(body(.rc_make_condition_celltype_metacells)), collapse = "\n"
  )
  native_text <- paste(
    deparse(body(.rc_native_supercell_membership)), collapse = "\n"
  )
  expect_match(body_text, "gamma = 30L", fixed = TRUE)
  expect_match(
    native_text, "SCimplify_by_graph_group_from_embedding", fixed = TRUE
  )
  expect_match(native_text, "cell.graph.group", fixed = TRUE)
  expect_match(native_text, "cell.split.condition", fixed = TRUE)
  expect_match(native_text, ".rc_scale_embedding_block_by_group", fixed = TRUE)
  expect_false(grepl("cell.annotation", native_text, fixed = TRUE))
  expect_false(grepl("supercell_stratum_col", body_text, fixed = TRUE))
  expect_false(grepl("stratum_col", body_text, fixed = TRUE))
})

test_that("stepwise meta-modules consume GRN and metacells", {
  f <- formals(rc_regcompass_step_meta_modules)
  expect_true(all(c("grn", "metacells", "gem", "outdir") %in% names(f)))
  expect_false("pando_args" %in% names(f))
})
