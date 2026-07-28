test_that("Pando extraction keeps only finite target-model R-squared values", {
  extraction <- paste(
    deparse(body(rc_extract_pando_tf_peak_gene)),
    collapse = "\n"
  )
  runner <- paste(
    deparse(body(.rc_run_condition_single_cell_grns_legacy)),
    collapse = "\n"
  )
  expect_match(extraction, "is.finite(rsq) & rsq >= min_model_rsq", fixed = TRUE)
  expect_match(extraction, "Pando target-model GOF", fixed = TRUE)
  expect_match(runner, "pando_evidence_filters", fixed = TRUE)
  expect_false(grepl("reliable_rsq", runner, fixed = TRUE))
  expect_false(exists(".rc_pando_rsq_is_reliable", inherits = TRUE))
})

test_that("Pando evidence-filter parameters are validated before inference", {
  expect_silent(.rc_validate_pando_evidence_filters(
    padj_threshold = 0.05,
    min_abs_estimate = 0,
    min_model_rsq = 0.1,
    require_padj = TRUE
  ))
  expect_error(
    .rc_validate_pando_evidence_filters(-0.1, 0, 0.1, TRUE),
    "padj_threshold"
  )
  expect_error(
    .rc_validate_pando_evidence_filters(1.1, 0, 0.1, TRUE),
    "padj_threshold"
  )
  expect_error(
    .rc_validate_pando_evidence_filters(0.05, -1, 0.1, TRUE),
    "min_abs_estimate"
  )
  expect_error(
    .rc_validate_pando_evidence_filters(0.05, 0, -0.1, TRUE),
    "min_model_rsq"
  )
  expect_error(
    .rc_validate_pando_evidence_filters(0.05, 0, 0.1, NA),
    "require_padj"
  )
  expect_error(
    .rc_run_condition_single_cell_grns(
      object = NULL,
      gem = list(model_info = list(species = "human")),
      outdir = tempfile(),
      genome = NULL,
      min_cells = 0,
      species = "human"
    ),
    "min_cells"
  )
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
