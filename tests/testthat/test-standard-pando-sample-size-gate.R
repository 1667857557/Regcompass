test_that("standard Pando uses fixed default correlation thresholds", {
  routed <- .rc_route_pando_infer_args(
    list(),
    condition_types = character(),
    standard_types = "T_cell"
  )

  expect_equal(routed$standard$tf_cor, 0.1)
  expect_equal(routed$standard$peak_cor, 0.05)
  expect_equal(routed$standard$adjust_method, "BH")
})

test_that("standard Pando preserves explicitly requested fixed thresholds", {
  routed <- .rc_route_pando_infer_args(
    list(tf_cor = 0.05, peak_cor = 0.02),
    condition_types = character(),
    standard_types = "T_cell"
  )

  expect_equal(routed$standard$tf_cor, 0.05)
  expect_equal(routed$standard$peak_cor, 0.02)
})

test_that("standard Pando single-process routing does not alter thresholds", {
  routed <- .rc_standard_pando_single_process_args(list(
    tf_cor = 0.05,
    peak_cor = 0.02,
    nthread = 8L
  ))

  expect_equal(routed$tf_cor, 0.05)
  expect_equal(routed$peak_cor, 0.02)
  expect_equal(routed$nthread, 1L)
})

test_that("condition routing remains independently configurable", {
  routed <- .rc_route_pando_infer_args(
    list(tf_cor = 0.2, peak_cor = 0.05),
    condition_types = "T_cell",
    standard_types = character()
  )

  expect_equal(routed$condition$tf_cor, 0.2)
  expect_equal(routed$condition$peak_cor, 0.05)
})
