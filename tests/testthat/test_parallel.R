test_that("rc_parallel_lapply supports forced sequential execution", {
  out <- rc_parallel_lapply(1:3, function(x) x + 1L, BPPARAM = FALSE)
  expect_equal(out, list(2L, 3L, 4L))
})

test_that("rc_default_bpparam can be disabled with worker option", {
  old <- getOption("RegCompassR.workers")
  on.exit(options(RegCompassR.workers = old), add = TRUE)
  options(RegCompassR.workers = 1)
  expect_null(rc_default_bpparam())
})
