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
  expect_equal(routed$standard$padj_threshold, 0.1)
  expect_equal(routed$condition$padj_threshold, 0.1)
  expect_identical(routed$condition$condition_ridge_control$cv_folds, 4L)
  resolved <- .rc_standard_pando_infer_args(routed$standard)
  expect_equal(resolved$padj_threshold, 0.1)
  expect_equal(attr(resolved, "regcompass_padj_threshold"), 0.1)
})

test_that("standard Pando defaults to the ridge K1 solver", {
  routed <- .rc_route_pando_infer_args(
    list(tf_cor = 0.05, peak_cor = 0.05),
    condition_types = character(),
    standard_types = "standard_celltype"
  )
  expect_identical(routed$standard$method, "ridge")
  expect_true(is.list(routed$standard$ridge_control))
  expect_identical(routed$standard$rank_action, "mark")
  expect_identical(routed$standard$min_residual_df, 1L)
  expect_equal(routed$standard$padj_threshold, 0.05)

  direct <- .rc_standard_pando_infer_args(list())
  expect_identical(direct$method, "ridge")
  expect_true(is.list(direct$ridge_control))
  expect_equal(direct$padj_threshold, .rc_standard_pando_padj_default)
  expect_equal(
    attr(direct, "regcompass_padj_threshold"),
    .rc_standard_pando_padj_default
  )
})

test_that("standard glm remains explicitly available with post-fit threshold", {
  routed <- .rc_route_pando_infer_args(
    list(method = "glm", tf_cor = 0.05, peak_cor = 0.05,
         padj_threshold = 0.01),
    condition_types = character(),
    standard_types = "standard_celltype"
  )
  expect_identical(routed$standard$method, "glm")
  expect_equal(routed$standard$padj_threshold, 0.01)
  expect_false("ridge_control" %in% names(routed$standard))

  resolved <- .rc_standard_pando_infer_args(routed$standard)
  expect_false("padj_threshold" %in% names(resolved))
  expect_equal(attr(resolved, "regcompass_padj_threshold"), 0.01)
})

test_that("standard padj threshold validates its range", {
  expect_error(
    .rc_standard_pando_infer_args(list(padj_threshold = 0)),
    "padj_threshold.*\\(0, 1\\)"
  )
  expect_error(
    .rc_standard_pando_infer_args(list(padj_threshold = 1)),
    "padj_threshold.*\\(0, 1\\)"
  )
})

test_that("ridge defaults live in canonical Stage 1 functions", {
  expect_false(exists(
    ".rc_route_pando_infer_args_before_default_ridge", inherits = TRUE
  ))
  expect_false(exists(
    ".rc_standard_pando_infer_args_before_default_ridge", inherits = TRUE
  ))
  route_text <- paste(deparse(body(.rc_route_pando_infer_args)), collapse = "\n")
  expect_match(route_text, 'standard_args$method %||% "ridge"', fixed = TRUE)
  standard_text <- paste(
    deparse(body(.rc_standard_pando_infer_args)), collapse = "\n"
  )
  expect_match(standard_text, 'args$method %||% "ridge"', fixed = TRUE)
})