test_that("canonical capacity uses one compiled path", {
  body_text <- paste(deparse(body(rc_reaction_capacity)), collapse = "\n")
  expect_true(grepl(".rc_gpr_capacity_cpp", body_text, fixed = TRUE))
  expect_false(grepl("rc_parallel_lapply", body_text, fixed = TRUE))
  source_text <- paste(
    readLines("R/gpr_capacity.R", warn = FALSE), collapse = "\n"
  )
  positions <- gregexpr(
    "rc_reaction_capacity <- function", source_text, fixed = TRUE
  )[[1L]]
  expect_identical(sum(positions > 0L), 1L)
})
