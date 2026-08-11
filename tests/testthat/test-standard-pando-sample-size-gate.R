test_that("standard Pando uses 0.05 default thresholds", {
  routed <- .rc_route_pando_infer_args(
    list(),
    condition_types = character(),
    standard_types = "T_cell"
  )

  expect_equal(routed$standard$tf_cor, 0.05)
  expect_equal(routed$standard$peak_cor, 0.05)
  expect_equal(routed$standard$padj_threshold, 0.05)
  expect_equal(routed$standard$adjust_method, "BH")
})

test_that("standard Pando preserves explicitly requested fixed thresholds", {
  routed <- .rc_route_pando_infer_args(
    list(tf_cor = 0.02, peak_cor = 0.03, padj_threshold = 0.01),
    condition_types = character(),
    standard_types = "T_cell"
  )

  expect_equal(routed$standard$tf_cor, 0.02)
  expect_equal(routed$standard$peak_cor, 0.03)
  expect_equal(routed$standard$padj_threshold, 0.01)
})

test_that("standard Pando single-process routing does not alter thresholds", {
  routed <- .rc_standard_pando_single_process_args(list(
    tf_cor = 0.02,
    peak_cor = 0.03,
    padj_threshold = 0.01,
    nthread = 8L
  ))

  expect_equal(routed$tf_cor, 0.02)
  expect_equal(routed$peak_cor, 0.03)
  expect_equal(routed$padj_threshold, 0.01)
  expect_equal(routed$nthread, 1L)
})

test_that("condition routing remains independently configurable", {
  routed <- .rc_route_pando_infer_args(
    list(tf_cor = 0.2, peak_cor = 0.04, padj_threshold = 0.02),
    condition_types = "T_cell",
    standard_types = character()
  )

  expect_equal(routed$condition$tf_cor, 0.2)
  expect_equal(routed$condition$peak_cor, 0.04)
  expect_equal(routed$condition$padj_threshold, 0.02)
})