test_that("single-condition standard Pando drops condition-only controls", {
  routed <- .rc_route_pando_infer_args(
    list(
      tf_cor = 0.1,
      peak_cor = 0,
      adjust_method = "BH",
      padj_threshold = 0.05,
      rank_action = "mark",
      min_residual_df = 1L
    ),
    condition_types = character(),
    standard_types = "T_cell"
  )

  expect_equal(routed$standard$tf_cor, 0.1)
  expect_equal(routed$standard$peak_cor, 0)
  expect_equal(routed$standard$adjust_method, "BH")
  expect_false(any(c(
    "padj_threshold", "rank_action", "min_residual_df"
  ) %in% names(routed$standard)))
  expect_setequal(
    routed$diagnostics$argument,
    c("padj_threshold", "rank_action", "min_residual_df")
  )
})

test_that("condition GRN drops standard-model controls", {
  routed <- .rc_route_pando_infer_args(
    list(
      tf_cor = 0.2,
      peak_cor = 0.05,
      adjust_method = "BH",
      padj_threshold = 0.05,
      rank_action = "mark",
      min_residual_df = 2L,
      method = "glmnet",
      alpha = 0.5,
      scale = TRUE
    ),
    condition_types = "Monocyte",
    standard_types = character()
  )

  expect_equal(routed$condition$tf_cor, 0.2)
  expect_equal(routed$condition$peak_cor, 0.05)
  expect_equal(routed$condition$min_residual_df, 2L)
  expect_false(any(c("method", "alpha", "scale") %in%
                     names(routed$condition)))
  expect_setequal(
    routed$diagnostics$argument,
    c("method", "alpha", "scale")
  )
})

test_that("mixed routing preserves controls for their own route", {
  routed <- .rc_route_pando_infer_args(
    list(
      tf_cor = 0.1,
      peak_cor = 0,
      adjust_method = "BH",
      padj_threshold = 0.05,
      rank_action = "mark",
      min_residual_df = 1L,
      method = "glm",
      scale = FALSE
    ),
    condition_types = "T_cell",
    standard_types = "HSPC"
  )

  expect_true(all(c("padj_threshold", "rank_action", "min_residual_df") %in%
                    names(routed$condition)))
  expect_true(all(c("method", "scale") %in% names(routed$standard)))
  expect_false("method" %in% names(routed$condition))
  expect_false("rank_action" %in% names(routed$standard))
})

test_that("unknown Pando arguments fail before model fitting", {
  expect_error(
    .rc_route_pando_infer_args(
      list(tf_cor = 0.1, misspelled_argument = TRUE),
      standard_types = "T_cell"
    ),
    "Unsupported `pando_infer_args`"
  )
})

test_that("canonical condition significance rule remains fixed", {
  expect_error(
    .rc_route_pando_infer_args(
      list(adjust_method = "BH", padj_threshold = 0.1),
      condition_types = "T_cell"
    ),
    "BH padj < 0.05"
  )
})
