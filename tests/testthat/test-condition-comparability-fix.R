test_that("condition comparison mask requires support in both conditions", {
  eligibility <- matrix(
    c(TRUE, TRUE, TRUE, FALSE),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("edge_both", "edge_reference_only"),
                    c("Control", "Drug"))
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

  observed <- .rc_condition_fit_comparison_mask(fit)

  expect_true(observed["edge_both", "Drug"])
  expect_false(observed["edge_reference_only", "Drug"])
  expect_true(all(observed[, "Control"] == eligibility[, "Control"]))
})

test_that("non-estimable reference contrasts cannot enter active effects", {
  eligibility <- matrix(
    c(TRUE, TRUE, TRUE, FALSE),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("edge_both", "edge_reference_only"),
                    c("Control", "Drug"))
  )
  comparison <- eligibility & matrix(
    eligibility[, "Control"],
    nrow = nrow(eligibility),
    ncol = ncol(eligibility),
    dimnames = dimnames(eligibility)
  )
  beta <- matrix(
    c(1, 2, 1, 0),
    nrow = 2L,
    byrow = TRUE,
    dimnames = dimnames(eligibility)
  )
  fit <- list(
    cell_type = "T",
    beta = beta,
    eligibility_mask = eligibility,
    comparison_mask = comparison,
    reference_condition = "Control"
  )
  condition_all <- data.frame(
    edge_id = rep(rownames(beta), times = 2L),
    condition = rep(colnames(beta), each = nrow(beta)),
    cell_type = "T",
    condition_estimate = c(beta[, "Control"], beta[, "Drug"]),
    condition_effect = c(0, 0, beta[, "Drug"] - beta[, "Control"]),
    estimate = c(beta[, "Control"], beta[, "Drug"]),
    rsq = 0.8,
    sample_blocked_oof_available = TRUE,
    eligible_in_condition = c(
      eligibility[, "Control"], eligibility[, "Drug"]
    ),
    stringsAsFactors = FALSE
  )
  effect_all <- condition_all
  effect_all$estimate <- effect_all$condition_effect
  extracted <- list(
    fit_contracts = list(fit),
    condition_all = condition_all,
    condition_active = condition_all,
    condition_effect_all = effect_all,
    condition_effect_active = effect_all,
    active_tol = 1e-8
  )

  observed <- .rc_apply_condition_comparison_semantics(
    extracted,
    condition_col = "condition",
    celltype_col = "cell_type",
    min_abs_estimate = 0,
    min_model_rsq = 0.1
  )

  drug_rows <- observed$condition_effect_all$condition == "Drug"
  comparable <- observed$condition_effect_all$comparable_to_reference[drug_rows]
  names(comparable) <- observed$condition_effect_all$edge_id[drug_rows]
  expect_true(comparable[["edge_both"]])
  expect_false(comparable[["edge_reference_only"]])
  expect_true("edge_both" %in% observed$condition_effect_active$edge_id)
  expect_false("edge_reference_only" %in%
                 observed$condition_effect_active$edge_id)

  single_sample <- extracted
  single_sample$condition_all$rsq <- NA_real_
  single_sample$condition_all$sample_blocked_oof_available <- FALSE
  single_sample$condition_effect_all$rsq <- NA_real_
  single_sample$condition_effect_all$sample_blocked_oof_available <- FALSE
  exploratory <- .rc_apply_condition_comparison_semantics(
    single_sample,
    condition_col = "condition",
    celltype_col = "cell_type",
    min_abs_estimate = 0,
    min_model_rsq = 0.1
  )
  expect_true("edge_both" %in% exploratory$condition_effect_active$edge_id)
  expect_match(
    exploratory$reliability_policy,
    "zero regulatory reliability",
    fixed = TRUE
  )
})

test_that("mouse analyses cannot silently use hg38 regulatory regions", {
  expect_error(
    .rc_default_pando_regions("mouse"),
    "hg38 conserved-element set is not valid for mouse input",
    fixed = TRUE
  )
})
