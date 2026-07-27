test_that("GREAT remains the canonical structural method", {
  expect_identical(
    .rc_validate_canonical_pando_design_args()$peak_to_gene_method,
    "GREAT"
  )
})
