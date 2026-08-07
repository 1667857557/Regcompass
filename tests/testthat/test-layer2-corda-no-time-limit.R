test_that("CORDA2 structural completion has no time limit", {
  expect_identical(
    RegCompassR:::.rc_layer2_completion_time_limit(list(), TRUE),
    Inf
  )

  expect_error(
    RegCompassR:::.rc_layer2_completion_time_limit(
      list(completion_time_limit = 3000), TRUE
    ),
    "CORDA2 reconstruction runs without a time limit"
  )

  implementation <- paste(
    deparse(body(RegCompassR:::.rc_layer2_prepare_completion)),
    collapse = "\n"
  )
  expect_match(
    implementation,
    "clean_params$completion_time_limit <- Inf",
    fixed = TRUE
  )
})

test_that("non-CORDA2 structural routes retain finite time-limit controls", {
  expect_equal(
    RegCompassR:::.rc_layer2_completion_time_limit(
      list(completion_time_limit = 1200), FALSE
    ),
    1200
  )
  expect_equal(
    RegCompassR:::.rc_layer2_completion_time_limit(list(), FALSE),
    300
  )
  expect_identical(
    RegCompassR:::.rc_layer2_completion_time_limit(
      list(completion_time_limit = Inf), FALSE
    ),
    Inf
  )
})
