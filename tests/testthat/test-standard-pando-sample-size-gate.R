test_that("standard Pando raises tf_cor for small cell groups", {
  gated <- .rc_standard_pando_sample_size_gate(
    list(tf_cor = 0.05, peak_cor = 0.05),
    n_cells = 300L
  )
  gate <- attr(gated, "sample_size_aware_tf_cor_gate", exact = TRUE)
  expected_t <- stats::qt(0.975, df = 298)
  expected_r <- expected_t / sqrt(expected_t^2 + 298)

  expect_equal(gated$tf_cor, expected_r, tolerance = 1e-12)
  expect_gt(gated$tf_cor, 0.05)
  expect_equal(gated$peak_cor, 0.05)
  expect_equal(gate$requested_tf_cor, 0.05)
  expect_equal(gate$sample_size_floor, expected_r, tolerance = 1e-12)
  expect_equal(gate$effective_tf_cor, expected_r, tolerance = 1e-12)
  expect_equal(gate$alpha, 0.05)
  expect_equal(gate$n_cells, 300L)
})

test_that("standard Pando keeps a stronger requested tf_cor floor", {
  gated <- .rc_standard_pando_sample_size_gate(
    list(tf_cor = 0.1, peak_cor = 0),
    n_cells = 5000L
  )
  gate <- attr(gated, "sample_size_aware_tf_cor_gate", exact = TRUE)

  expect_equal(gated$tf_cor, 0.1)
  expect_lt(gate$sample_size_floor, 0.1)
  expect_equal(gate$effective_tf_cor, 0.1)
  expect_equal(gated$peak_cor, 0)
})

test_that("standard Pando uses the original tf_cor default as effect-size floor", {
  gated <- .rc_standard_pando_sample_size_gate(list(), n_cells = 5000L)
  gate <- attr(gated, "sample_size_aware_tf_cor_gate", exact = TRUE)

  expect_equal(gate$requested_tf_cor, 0.1)
  expect_equal(gated$tf_cor, 0.1)
})

test_that("sample-size-aware gate rejects invalid inputs", {
  expect_error(
    .rc_standard_pando_sample_size_gate(list(tf_cor = 0.1), n_cells = 3L),
    "n_cells"
  )
  expect_error(
    .rc_standard_pando_sample_size_gate(list(tf_cor = 1), n_cells = 300L),
    "tf_cor"
  )
})

test_that("condition routing remains configurable and is not sample-size gated", {
  routed <- .rc_route_pando_infer_args(
    list(tf_cor = 0.2, peak_cor = 0.05),
    condition_types = "T_cell",
    standard_types = character()
  )

  expect_equal(routed$condition$tf_cor, 0.2)
  expect_equal(routed$condition$peak_cor, 0.05)
})
