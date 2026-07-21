test_that("missing Pando rsq is never treated as reliable", {
  expect_identical(
    .rc_pando_rsq_is_reliable(
      c(NA_real_, NaN, Inf, -Inf, 0, 0.1, 0.2),
      min_model_rsq = 0.1
    ),
    c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE)
  )
  expect_false(.rc_pando_rsq_is_reliable(NA_real_, min_model_rsq = 0))
})

test_that("zero regulatory reliability falls back to RNA support", {
  rna <- matrix(
    c(0, 0.2, 0.8, 1),
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("m1", "m2"))
  )
  modifier <- matrix(0, nrow = 2, ncol = 2, dimnames = dimnames(rna))
  integrated <- .rc_integrate_regulatory_support_v170(
    rna_support = rna,
    regulatory_modifier = modifier,
    alpha = 1
  )
  expect_equal(as.numeric(integrated), as.numeric(rna), tolerance = 0)
  expect_identical(dimnames(integrated), dimnames(rna))
  expect_match(
    attr(integrated, "integration_formula"),
    "C_multiome",
    fixed = TRUE
  )
  expect_match(
    attr(integrated, "score_semantics"),
    "zero-preserving",
    fixed = TRUE
  )
})
