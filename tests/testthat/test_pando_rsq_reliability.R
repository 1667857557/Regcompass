test_that("Pando extraction keeps only finite target-model R-squared values", {
  extraction <- paste(
    deparse(body(rc_extract_pando_tf_peak_gene)),
    collapse = "\n"
  )
  runner <- paste(
    deparse(body(.rc_run_condition_single_cell_grns)),
    collapse = "\n"
  )
  expect_match(extraction, "is.finite(rsq) & rsq >= min_model_rsq", fixed = TRUE)
  expect_match(extraction, "Pando target-model GOF", fixed = TRUE)
  expect_false(grepl("reliable_rsq", runner, fixed = TRUE))
  expect_false(exists(".rc_pando_rsq_is_reliable", inherits = TRUE))
})

test_that("zero regulatory modifier falls back to RNA support", {
  rna <- matrix(
    c(0, 0.2, 0.8, 1),
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("m1", "m2"))
  )
  modifier <- matrix(0, nrow = 2, ncol = 2, dimnames = dimnames(rna))
  integrated <- .rc_integrate_regulatory_support(
    rna_support = rna,
    regulatory_modifier = modifier,
    alpha = 1
  )
  expect_equal(as.numeric(integrated), as.numeric(rna), tolerance = 0)
  expect_identical(dimnames(integrated), dimnames(rna))
  expect_match(attr(integrated, "integration_formula"), "C_multiome", fixed = TRUE)
  expect_match(attr(integrated, "score_semantics"), "zero-preserving", fixed = TRUE)
})
