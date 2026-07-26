test_that("training-fold condition centers do not use validation cells", {
  condition <- c(rep("A", 4), rep("B", 4))
  x <- cbind(
    x1 = c(1, 2, 3, 1000, 5, 6, 7, -1000),
    x2 = c(2, 4, 6, 2000, 1, 3, 5, -2000)
  )
  y <- c(1, 2, 3, 1000, 5, 6, 7, -1000)
  training <- c(1:3, 5:7)
  validation <- c(4, 8)

  x_centers <- RegCompassR:::.rc_condition_matrix_centers(
    x[training, , drop = FALSE], condition[training]
  )
  y_centers <- RegCompassR:::.rc_condition_vector_centers(
    y[training], condition[training]
  )

  expect_equal(x_centers["A", ], c(2, 4))
  expect_equal(x_centers["B", ], c(6, 3))
  expect_equal(y_centers[c("A", "B")], c(A = 2, B = 6))

  x_validation <- RegCompassR:::.rc_apply_matrix_centers(
    x[validation, , drop = FALSE],
    condition[validation],
    x_centers
  )
  y_validation <- RegCompassR:::.rc_apply_vector_centers(
    y[validation], condition[validation], y_centers
  )
  expect_equal(x_validation[1, ], c(998, 1996))
  expect_equal(x_validation[2, ], c(-1006, -2003))
  expect_equal(y_validation, c(998, -1006))
})

test_that("candidate outcome screening is disabled in canonical mode", {
  expect_error(
    RegCompassR:::.rc_validate_multitask_grn_args(list(
      candidate_screen_threshold = 0.01
    )),
    "must be 0"
  )
  args <- RegCompassR:::.rc_validate_multitask_grn_args(list(
    candidate_screen_threshold = 0
  ))
  expect_identical(args$candidate_screen_threshold, 0)
})

test_that("manual CV records fold-specific preprocessing", {
  expect_true(is.function(RegCompassR:::.rc_multitask_manual_cv))
  expect_true(is.function(RegCompassR:::.rc_condition_matrix_centers))
  expect_true(is.function(RegCompassR:::.rc_equal_condition_edge_scale))
  text <- paste(
    deparse(body(RegCompassR:::.rc_multitask_manual_cv)),
    collapse = "\n"
  )
  expect_match(text, "training_condition", fixed = TRUE)
  expect_match(text, "x_centers", fixed = TRUE)
  expect_match(text, "train_scale", fixed = TRUE)
  expect_false(grepl("cv.glmnet", text, fixed = TRUE))
})

test_that("edge ordering is fixed before predictor matrices are constructed", {
  text <- paste(
    deparse(body(RegCompassR:::.rc_fit_multitask_target)),
    collapse = "\n"
  )
  expect_match(text, "order(as.character(edges$edge_id))", fixed = TRUE)
})
