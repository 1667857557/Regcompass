test_that("cached union GEM reuse requires checksum and medium identity", {
  S <- matrix(
    0,
    nrow = 1,
    ncol = 1,
    dimnames = list("M1", "R1")
  )
  model <- rc_make_gem(S, lb = c(R1 = 0), ub = c(R1 = 1))
  model$is_union_gem <- TRUE
  model$union_gem_medium_scenario <- "medium_a"

  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file), add = TRUE)
  saveRDS(model, file)
  checksum <- unname(tools::md5sum(file))

  expect_silent(.rc_read_cached_union_gem(
    file = file,
    medium_scenario = "medium_a",
    expected_checksum = checksum
  ))
  expect_error(
    .rc_read_cached_union_gem(
      file = file,
      medium_scenario = "medium_a",
      expected_checksum = ""
    ),
    "checksum is missing"
  )
  expect_error(
    .rc_read_cached_union_gem(
      file = file,
      medium_scenario = "medium_a",
      expected_checksum = "incorrect"
    ),
    "checksum check"
  )
  expect_error(
    .rc_read_cached_union_gem(
      file = file,
      medium_scenario = "medium_b",
      expected_checksum = checksum
    ),
    "original final medium-specific union GEM"
  )
})

test_that("plain GEM files cannot be reused as final union GEMs", {
  S <- matrix(
    0,
    nrow = 1,
    ncol = 1,
    dimnames = list("M1", "R1")
  )
  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file), add = TRUE)
  saveRDS(rc_make_gem(S, lb = c(R1 = 0), ub = c(R1 = 1)), file)

  expect_error(
    .rc_read_cached_union_gem(
      file = file,
      medium_scenario = "medium_a",
      expected_checksum = unname(tools::md5sum(file))
    ),
    "original final medium-specific union GEM"
  )
})
