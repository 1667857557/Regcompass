test_that("runtime RDS storage is gzip-compressed", {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path, force = TRUE), add = TRUE)
  value <- list(integer = 1:5, text = "RegCompass")

  RegCompassR:::saveRDS(value, path)

  expect_true(file.exists(path))
  expect_identical(
    readBin(path, what = "raw", n = 2L),
    as.raw(c(0x1f, 0x8b))
  )
  expect_identical(readRDS(path), value)
})

test_that("duplicate Stage 1 GRN RDS is not written", {
  directory <- tempfile("regcompass-rds-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE, force = TRUE), add = TRUE)
  path <- file.path(directory, "single_cell_grn.rds")

  RegCompassR:::saveRDS(list(stage = 1L), path)

  expect_false(file.exists(path))
})
