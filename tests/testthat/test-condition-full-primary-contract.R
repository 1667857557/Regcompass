test_that("Layer 1 uses the padj-filtered fixed-dictionary projection", {
  helper <- paste(
    deparse(body(.rc_condition_pando_projection)), collapse = "\n"
  )
  expect_match(helper, "Pando::project_condition_grn_cells", fixed = TRUE)
  expect_match(helper, "significant_only = TRUE", fixed = TRUE)
  expect_match(
    helper,
    "Pando::aggregate_condition_grn_projection",
    fixed = TRUE
  )
  expect_match(
    helper,
    "paired_cell_full_fit_fixed_dictionary_glm_padj_filtered",
    fixed = TRUE
  )
  expect_false(grepl("project_condition_grn_primary_cells", helper, fixed = TRUE))
  expect_false(grepl("pairwise_common", helper, fixed = TRUE))
  expect_false(grepl("global_common", helper, fixed = TRUE))
})

test_that("condition compatibility decomposition is primary plus zero", {
  helper <- paste(
    deparse(body(.rc_condition_pando_projection)), collapse = "\n"
  )
  expect_match(helper, "common = primary", fixed = TRUE)
  expect_match(helper, "condition_unique = condition_unique", fixed = TRUE)
  expect_match(helper, "primary * 0", fixed = TRUE)
})

test_that("Layer 1 validator accepts new provenance and aliases", {
  validator <- paste(
    deparse(body(.rc_validate_layer1_stage)), collapse = "\n"
  )
  expect_match(
    validator,
    "paired_cell_full_fit_fixed_dictionary_glm_padj_filtered",
    fixed = TRUE
  )
  expect_match(
    validator,
    "padj_filtered_fixed_dictionary_condition_glm",
    fixed = TRUE
  )
  expect_match(
    validator,
    "coefficient_NA_and_zero_realized_penalty_contribution",
    fixed = TRUE
  )
  expect_false(grepl("outer_condition_stratified_cell_oof", validator,
                     fixed = TRUE))
  expect_false(grepl("structural_zero_by_condition", validator,
                     fixed = TRUE))
})

test_that("Layer 2 keeps shared model comparability", {
  text <- paste(deparse(body(rc_regcompass_step_layer2)), collapse = "\n")
  expect_match(text, "penalty_condition_full_oof", fixed = TRUE)
  expect_match(text, "penalty_common_oof", fixed = TRUE)
  expect_match(text, "penalty_condition_unique_increment", fixed = TRUE)
  expect_match(text, "penalty_rna_only", fixed = TRUE)
  removed <- c(
    "penalty_depth_matched_rna",
    "penalty_common_depth_interval_rna",
    "penalty_alpha_sensitivity",
    "zero_support_sensitive",
    "link_saturation_sensitive"
  )
  expect_false(any(vapply(removed, grepl, logical(1), x = text, fixed = TRUE)))
})
