test_that("standard Pando candidate intervals map back to ATAC peaks", {
  map <- RegCompassR:::.rc_map_pando_regions_to_atac

  expect_identical(
    map(
      c("chr1:120-180", "chr2-50-75", "chr3:1-10"),
      c("chr1-100-200", "chr2:25-100")
    ),
    c(1L, 2L, NA_integer_)
  )
})

test_that("standard Pando peak mapping is deterministic for overlaps", {
  map <- RegCompassR:::.rc_map_pando_regions_to_atac

  # The first candidate shares more bases with the second peak. The second
  # candidate has an equal overlap and midpoint distance, so assay order wins.
  expect_identical(
    map(
      c("chr1-175-240", "chr2-150-150"),
      c("chr1-100-200", "chr1-180-260", "chr2-99-199", "chr2-101-201")
    ),
    c(2L, 3L)
  )
})

test_that("standard Pando peak mapping rejects normalized duplicates", {
  expect_error(
    RegCompassR:::.rc_map_pando_regions_to_atac(
      "chr1-120-180", c("chr1:100-200", "chr1-100-200")
    ),
    "duplicated after Pando region normalization"
  )
})
