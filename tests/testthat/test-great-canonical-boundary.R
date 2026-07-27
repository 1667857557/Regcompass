test_that("canonical runner validates GREAT before calling Pando", {
  defaults <- .rc_validate_canonical_pando_design_args(list())
  expect_identical(defaults$peak_to_gene_method, "GREAT")

  signac <- .rc_validate_canonical_pando_design_args(list(
    peak_to_gene_method = "Signac"
  ))
  expect_identical(signac$peak_to_gene_method, "Signac")
})
