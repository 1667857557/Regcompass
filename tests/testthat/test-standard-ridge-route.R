test_that("standard ridge routes its controls", {
  routed <- .rc_route_pando_infer_args(
    list(
      method = "ridge",
      tf_cor = 0.1,
      peak_cor = 0.05,
      padj_threshold = 0.1,
      rank_action = "mark",
      min_residual_df = 2L,
      ridge_control = list(cv_folds = 3L),
      condition_ridge_control = list(cv_folds = 4L)
    ),
    condition_types = "condition_celltype",
    standard_types = "standard_celltype"
  )
  expect_identical(routed$standard$method, "ridge")
  expect_identical(routed$standard$ridge_control$cv_folds, 3L)
  expect_identical(routed$standard$rank_action, "mark")
  expect_identical(routed$standard$min_residual_df, 2L)
  expect_false("padj_threshold" %in% names(routed$standard))
  expect_equal(routed$condition$padj_threshold, 0.1)
  expect_identical(routed$condition$condition_ridge_control$cv_folds, 4L)
  resolved <- .rc_standard_pando_infer_args(routed$standard)
  expect_equal(resolved$padj_threshold, .rc_standard_pando_padj_fixed)
})

test_that("standard glm remains available", {
  routed <- .rc_route_pando_infer_args(
    list(method = "glm", tf_cor = 0.1, peak_cor = 0.05),
    condition_types = character(),
    standard_types = "standard_celltype"
  )
  expect_identical(routed$standard$method, "glm")
  expect_false("ridge_control" %in% names(routed$standard))
})

test_that("standard ridge is implemented in canonical source functions", {
  expect_false(exists(
    ".rc_pando_infer_arg_catalog_standard_ridge", inherits = TRUE
  ))
  expect_false(exists(
    ".rc_route_pando_infer_args_standard_ridge", inherits = TRUE
  ))
  expect_false(exists(
    ".rc_run_standard_pando_celltype_job_ridge", inherits = TRUE
  ))
  expect_false(exists(
    ".rc_fit_pando_by_celltype_route_ridge", inherits = TRUE
  ))
  route_text <- paste(deparse(body(.rc_route_pando_infer_args)), collapse = "\n")
  expect_match(route_text, "ridge_control", fixed = TRUE)
  standard_text <- paste(
    deparse(body(.rc_standard_pando_infer_args)), collapse = "\n"
  )
  expect_match(standard_text, "is_ridge", fixed = TRUE)
})
