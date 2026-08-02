common_dictionary_fit_fixture <- function() {
  edge <- data.frame(
    edge_id = c("G||TF1||P1", "G||TF2||P2"),
    target = "G",
    tf = c("TF1", "TF2"),
    region = c("P1", "P2"),
    atac_feature_id = c("A1", "A2"),
    candidate_index = 1:2,
    stringsAsFactors = FALSE
  )
  coefficient <- do.call(rbind, lapply(c("A", "B"), function(condition) {
    data.frame(
      edge_id = edge$edge_id,
      target = edge$target,
      tf = edge$tf,
      region = edge$region,
      condition = condition,
      estimate = c(1, -0.5),
      std_err = c(0.1, 0.2),
      statistic = c(10, -2.5),
      pval = c(1e-6, 0.02),
      padj = c(2e-6, 0.08),
      significant = c(TRUE, FALSE),
      penalty_effect = c(1, 0),
      estimable = TRUE,
      zero_variance = FALSE,
      aliased = FALSE,
      stringsAsFactors = FALSE
    )
  }))
  structure(list(
    schema_version = "pando_condition_grn_common_dictionary_v1",
    fit_engine = "two_stage_exact_edge_union_fixed_dictionary_glm",
    coefficient_scale = "shared_preprocessed_input_units_unscaled",
    inference_scope = "conditional_on_selected_edge_dictionary",
    cell_type = "T_cell",
    condition_levels = c("A", "B"),
    condition_col = "condition",
    cell_type_col = "cell_type",
    condition_cell_ids = list(A = c("a1", "a2"), B = c("b1", "b2")),
    edge_dictionary = edge,
    coefficients = coefficient,
    fit = data.frame(
      target = rep("G", 2), condition = c("A", "B"),
      rsq = c(0.8, 0.7), fit_status = "ok",
      stringsAsFactors = FALSE
    ),
    network_names = c(A = "net_A", B = "net_B"),
    padj_threshold = 0.05,
    adjust_method = "BH",
    scale = FALSE,
    interaction = ":",
    projection_effect_column = "penalty_effect",
    projection_policy = "padj_significant_effects_only"
  ), class = c("ConditionGRNFit", "list"))
}

test_that("RegCompass accepts only the common-dictionary Pando schema", {
  canonical <- common_dictionary_fit_fixture()
  expect_invisible(
    RegCompassR:::.rc_require_pando_condition_grn_fit(canonical)
  )
  legacy <- canonical
  legacy$schema_version <- "pando_condition_grn_fit"
  expect_error(
    RegCompassR:::.rc_require_pando_condition_grn_fit(legacy),
    "common-dictionary"
  )
})

test_that("penalty effect is exactly the BH-significant coefficient", {
  invalid <- common_dictionary_fit_fixture()
  invalid$coefficients$penalty_effect[
    !invalid$coefficients$significant
  ] <- 0.2
  expect_error(
    RegCompassR:::.rc_require_pando_condition_grn_fit(invalid),
    "invalid"
  )

  invalid <- common_dictionary_fit_fixture()
  invalid$coefficients$penalty_effect[
    invalid$coefficients$significant
  ] <- 0
  expect_error(
    RegCompassR:::.rc_require_pando_condition_grn_fit(invalid),
    "invalid"
  )
})

test_that("retired condition engine is absent from package metadata", {
  root <- testthat::test_path("..", "..")
  description <- read.dcf(file.path(root, "DESCRIPTION"))
  expect_identical(
    unname(description[1L, "Config/RegCompass/PandoConditionSchema"]),
    "pando_condition_grn_common_dictionary_v1"
  )
  expect_identical(
    unname(description[1L, "Config/RegCompass/PandoConditionEffectFilter"]),
    "BH-adjusted-p-below-0.05"
  )
  expect_false(file.exists(file.path(root, "R", "condition_grn_runtime.R")))
  expect_false(file.exists(file.path(root, "R", "condition_grn_runtime_guard.R")))
})
