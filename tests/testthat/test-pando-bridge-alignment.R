test_that("Pando comparison support must be explicit and reference aligned", {
  eligibility <- matrix(
    c(TRUE, TRUE, TRUE, FALSE),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(
      c("edge_both", "edge_reference_only"),
      c("Control", "Drug")
    )
  )
  comparison <- eligibility & matrix(
    eligibility[, "Control"],
    nrow = nrow(eligibility),
    ncol = ncol(eligibility),
    dimnames = dimnames(eligibility)
  )
  fit <- list(
    beta = matrix(
      c(1, 2, 1, 0),
      nrow = 2L,
      byrow = TRUE,
      dimnames = dimnames(eligibility)
    ),
    eligibility_mask = eligibility,
    comparison_mask = comparison,
    reference_condition = "Control"
  )

  expect_identical(.rc_condition_fit_comparison_mask(fit), comparison)
  fit$comparison_mask <- NULL
  expect_error(
    .rc_condition_fit_comparison_mask(fit),
    "Pando_regcompass >= 1.2.1",
    fixed = TRUE
  )
})

test_that("RegCompass protects the Pando bridge contract", {
  expect_error(
    .rc_validate_pando_bridge_args(
      pando_infer_args = list(genes = "HK2")
    ),
    "RegCompass-managed fields: genes",
    fixed = TRUE
  )
  expect_error(
    .rc_validate_pando_bridge_args(
      pando_initiate_args = list(object = "wrong")
    ),
    "RegCompass-managed fields: object",
    fixed = TRUE
  )
  expect_error(
    .rc_validate_pando_bridge_args(
      pando_motif_args = list(genome = "wrong")
    ),
    "RegCompass-managed fields: genome",
    fixed = TRUE
  )
  expect_error(
    .rc_validate_pando_bridge_args(
      pando_infer_args = list(aggregate_rna_col = "metacell")
    ),
    "Stage 2 owns metacell aggregation",
    fixed = TRUE
  )
  expect_silent(.rc_validate_pando_bridge_args(
    pando_infer_args = list(
      candidate_screen = "motif_domain",
      reference_condition = "Control"
    )
  ))
})

test_that("canonical Stage 1 defaults and parallel routing match Pando", {
  infer_defaults <- eval(
    formals(.rc_run_condition_single_cell_grns)$pando_infer_args
  )
  expect_identical(infer_defaults$method, "shared_design_independent")
  expect_identical(infer_defaults$candidate_screen, "motif_domain")
  expect_identical(infer_defaults$condition_mix, 1)
  expect_identical(infer_defaults$condition_weight, "equal")
  expect_true(infer_defaults$scale)

  step_body <- paste(deparse(body(rc_regcompass_step_grn)), collapse = "\n")
  expect_match(step_body, "use_bpparam", fixed = TRUE)
  expect_match(step_body, "use_pando_native", fixed = TRUE)
  expect_match(step_body, "identical(BPPARAM, TRUE)", fixed = TRUE)
  expect_match(step_body, "pando_parallel", fixed = TRUE)
})

test_that("Pando bridge provenance records the active policy", {
  run_body <- paste(
    deparse(body(.rc_run_condition_single_cell_grns)),
    collapse = "\n"
  )
  expect_match(run_body, "pando_candidate_screen", fixed = TRUE)
  expect_match(run_body, "comparison_support", fixed = TRUE)
  expect_match(run_body, "single_cell_grn.rds", fixed = TRUE)
})
