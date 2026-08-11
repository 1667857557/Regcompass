.strict_fit_fixture <- function(padj_threshold = 0.05) {
  edge <- data.frame(
    edge_id = c("G||TF1||P1", "G||TF2||P2"),
    target = "G",
    tf = c("TF1", "TF2"),
    region = c("P1", "P2"),
    atac_feature_id = c("A1", "A2"),
    candidate_index = 1:2,
    screening_min_padj = c(0.002, 0.04),
    screening_n_significant_conditions = c(2L, 1L),
    screening_significant_conditions = c("A;B", "A"),
    stringsAsFactors = FALSE
  )
  attr(edge, "preprocessing_provenance_verified") <- TRUE

  pval <- c(0.001, 0.04, 0.02, 0.8)
  condition <- rep(c("A", "B"), each = 2L)
  padj <- unlist(lapply(split(pval, condition), stats::p.adjust,
                        method = "BH"), use.names = FALSE)
  estimate <- c(1, -0.2, -0.7, 0.1)
  shared <- c(0.15, -0.05, 0.15, -0.05)
  significant <- padj < padj_threshold
  coefficient <- data.frame(
    edge_id = rep(edge$edge_id, 2L),
    target = "G",
    tf = rep(edge$tf, 2L),
    region = rep(edge$region, 2L),
    condition = condition,
    estimate = estimate,
    shared_estimate = shared,
    condition_deviation = estimate - shared,
    std_err = 0.1,
    statistic = estimate / 0.1,
    pval = pval,
    padj = padj,
    significant = significant,
    penalty_effect = ifelse(significant, estimate, 0),
    estimable = TRUE,
    zero_variance = FALSE,
    aliased = FALSE,
    direction = ifelse(estimate > 0, "positive", "negative"),
    stringsAsFactors = FALSE
  )

  contrast_pval <- c(0.01, 0.5)
  contrast_padj <- stats::p.adjust(contrast_pval, method = "BH")
  contrasts <- data.frame(
    edge_id = edge$edge_id,
    target = "G",
    tf = edge$tf,
    region = edge$region,
    atac_feature_id = edge$atac_feature_id,
    condition_a = "A",
    condition_b = "B",
    contrast = "A-B",
    estimate_a = estimate[1:2],
    estimate_b = estimate[3:4],
    contrast_estimate = estimate[1:2] - estimate[3:4],
    contrast_se = c(0.2, 0.2),
    contrast_statistic = c(8.5, -1.5),
    contrast_pval = contrast_pval,
    contrast_padj = contrast_padj,
    contrast_estimable = TRUE,
    contrast_significant = contrast_padj < padj_threshold,
    stringsAsFactors = FALSE
  )

  screening <- data.frame(
    edge_id = c(edge$edge_id, "G||TF3||P3"),
    screening_min_padj = c(0.002, 0.04, 0.4),
    screening_n_significant_conditions = c(2L, 1L, 0L),
    screening_significant_conditions = c("A;B", "A", ""),
    stringsAsFactors = FALSE
  )

  structure(list(
    schema_version = "pando_condition_grn_common_dictionary_v1",
    model_schema = "pando_condition_grn_multitask_ridge_v2",
    fit_engine = "two_stage_exact_edge_union_multitask_ridge",
    coefficient_scale = "raw_tf_atac_interaction_units",
    internal_predictor_scale = "equal_condition_within_condition_rms",
    inference_scope =
      "approximate_ridge_wald_diagnostic_conditional_on_dictionary_cv_lambda_and_fusion",
    cell_type = "T_cell",
    condition_levels = c("A", "B"),
    condition_col = "condition",
    cell_type_col = "cell_type",
    condition_cell_ids = list(A = c("a1", "a2"), B = c("b1", "b2")),
    edge_dictionary = edge,
    coefficients = coefficient,
    contrasts = contrasts,
    fit = data.frame(
      target = rep("G", 2), condition = c("A", "B"),
      rsq = c(0.8, 0.7), rsq_oof = c(0.8, 0.7),
      fit_status = "ok", lambda = 0.1,
      predictor_scale_reference = "equal_condition_within_condition_rms",
      stringsAsFactors = FALSE
    ),
    network_names = c(A = "net_A", B = "net_B"),
    padj_threshold = padj_threshold,
    adjust_method = "BH",
    scale = FALSE,
    interaction = ":",
    projection_effect_column = "penalty_effect",
    projection_policy = "padj_significant_ridge_effects",
    fit_dictionary_policy =
      "preliminary_joint_ridge_bh_significant_union_then_joint_refit",
    candidate_edge_count = 3L,
    fit_dictionary_edge_count = 2L,
    dictionary_screening_threshold = padj_threshold,
    dictionary_screening_summary = screening,
    rna_layer = "data",
    peak_layer = "data",
    peak_value_type = "normalized",
    preprocessing_fingerprint = "fixture-preprocessing",
    target_genes = "G"
  ), class = c("ConditionGRNFit", "list"))
}

test_that("significant-union dictionary is complete in every condition", {
  fit <- .strict_fit_fixture()
  expect_invisible(RegCompassR:::.rc_require_pando_condition_grn_fit(fit))
  incomplete <- fit
  incomplete$coefficients <- incomplete$coefficients[-1L, , drop = FALSE]
  expect_error(
    RegCompassR:::.rc_require_pando_condition_grn_fit(incomplete),
    "every screened fit-dictionary edge"
  )
})

test_that("every fit-dictionary edge must have preliminary BH support", {
  fit <- .strict_fit_fixture()
  bad <- fit
  bad$edge_dictionary$screening_n_significant_conditions[[2L]] <- 0L
  bad$edge_dictionary$screening_min_padj[[2L]] <- 0.4
  expect_error(
    RegCompassR:::.rc_require_pando_condition_grn_fit(bad),
    "BH-significant in at least one condition"
  )
})

test_that("final condition penalty uses configured BH threshold", {
  fit <- .strict_fit_fixture()
  expect_invisible(RegCompassR:::.rc_require_pando_condition_grn_fit(fit))

  nonsignificant <- which(!fit$coefficients$significant)[[1L]]
  expect_equal(fit$coefficients$penalty_effect[[nonsignificant]], 0)

  wrong_effect <- fit
  wrong_effect$coefficients$penalty_effect[[nonsignificant]] <-
    wrong_effect$coefficients$estimate[[nonsignificant]]
  expect_error(
    RegCompassR:::.rc_require_pando_condition_grn_fit(wrong_effect),
    "BH-significant finite ridge coefficient"
  )

  wrong_flag <- fit
  wrong_flag$coefficients$significant[[1L]] <- FALSE
  expect_error(
    RegCompassR:::.rc_require_pando_condition_grn_fit(wrong_flag),
    "configured BH threshold"
  )
})

test_that("padj threshold is configurable above 0.1 with no artificial cap", {
  fit <- .strict_fit_fixture(padj_threshold = 0.2)
  expect_invisible(RegCompassR:::.rc_require_pando_condition_grn_fit(fit))

  invalid <- fit
  invalid$padj_threshold <- 1
  invalid$dictionary_screening_threshold <- 1
  expect_error(
    RegCompassR:::.rc_require_pando_condition_grn_fit(invalid),
    "incomplete"
  )
})

test_that("pairwise differential GRN contrasts remain separate inference", {
  fit <- .strict_fit_fixture()
  expect_invisible(RegCompassR:::.rc_require_pando_condition_grn_fit(fit))
  expect_equal(
    fit$contrasts$contrast_estimate,
    fit$contrasts$estimate_a - fit$contrasts$estimate_b
  )

  wrong <- fit
  wrong$contrasts$contrast_padj[[1L]] <- 0.99
  expect_error(
    RegCompassR:::.rc_require_pando_condition_grn_fit(wrong),
    "contrast padj"
  )
})

test_that("Layer 1 keeps only final BH-significant ridge effects", {
  fit <- .strict_fit_fixture()
  gated <- RegCompassR:::.rc_apply_condition_penalty_gate(fit)
  threshold <- fit$padj_threshold
  expected <- fit$coefficients$estimable &
    is.finite(fit$coefficients$padj) &
    fit$coefficients$padj < threshold

  expect_identical(gated$coefficients$significant, expected)
  expect_equal(
    gated$coefficients$penalty_effect,
    ifelse(expected, gated$coefficients$estimate, 0)
  )
  expect_true(all(gated$coefficients$fit_status == "ok"))
  expect_identical(
    gated$regcompass_penalty_filter,
    "estimable & finite estimate & fit_status == 'ok' & BH padj < 0.05"
  )
  expect_identical(
    gated$regcompass_significance_role,
    "condition_edge_inclusion_gate"
  )
  expect_invisible(
    RegCompassR:::.rc_require_layer1_condition_grn_fit(gated)
  )
})

test_that("non-estimable condition edges contribute zero", {
  fit <- .strict_fit_fixture()
  fit$coefficients$estimable[[2L]] <- FALSE
  fit$coefficients$penalty_effect[[2L]] <- 0
  fit$coefficients$significant[[2L]] <- FALSE
  fit$coefficients$padj[[2L]] <- NA_real_
  fit$coefficients$pval[[2L]] <- NA_real_

  gated <- RegCompassR:::.rc_apply_condition_penalty_gate(fit)
  expect_equal(gated$coefficients$penalty_effect[[2L]], 0)
  expect_false(gated$coefficients$significant[[2L]])
  expect_invisible(
    RegCompassR:::.rc_require_layer1_condition_grn_fit(gated)
  )
})

test_that("target reliability requires at least one significant edge", {
  fit <- .strict_fit_fixture()
  reliability <- RegCompassR:::.rc_condition_target_reliability(fit)
  expect_true(all(reliability$n_significant_edges > 0L))
  expect_true(all(is.finite(reliability$reliability)))

  none <- fit
  none$coefficients$significant[none$coefficients$condition == "B"] <- FALSE
  none$coefficients$padj[none$coefficients$condition == "B"] <- 1
  none$coefficients$penalty_effect[none$coefficients$condition == "B"] <- 0
  reliability_none <- RegCompassR:::.rc_condition_target_reliability(none)
  expect_true(is.na(reliability_none$reliability[
    reliability_none$condition == "B"
  ]))
})

test_that("Layer 1 keeps product-of-means projection with BH gate", {
  selector <- paste(
    deparse(body(RegCompassR:::.rc_condition_pando_object_for_fit)),
    collapse = "\n"
  )
  projection <- paste(
    deparse(body(RegCompassR:::.rc_condition_pando_projection)),
    collapse = "\n"
  )
  expect_match(
    selector, ".rc_require_layer1_condition_grn_fit", fixed = TRUE
  )
  expect_match(
    projection, ".rc_require_layer1_condition_grn_fit", fixed = TRUE
  )
  expect_match(
    projection, ".rc_condition_penalty_gate", fixed = TRUE
  )
  expect_match(
    projection,
    "beta_times_group_mean_tf_times_group_mean_atac",
    fixed = TRUE
  )
  expect_false(grepl("mean(TF * ATAC)", projection, fixed = TRUE))
})

test_that("condition contract is implemented once in functional source files", {
  root <- testthat::test_path("..", "..")
  r_files <- list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)
  source_text <- unlist(lapply(r_files, readLines, warn = FALSE), use.names = FALSE)
  count_definition <- function(name) {
    pattern <- paste0("^", gsub("\\.", "\\\\.", name),
                      "[[:space:]]*<-[[:space:]]*function")
    sum(grepl(pattern, source_text))
  }
  expect_identical(count_definition(".rc_require_pando_condition_grn_fit"), 1L)
  expect_identical(count_definition(".rc_require_layer1_condition_grn_fit"), 1L)
  expect_identical(count_definition(".rc_extract_condition_grn_contract"), 1L)
  expect_identical(count_definition(".rc_condition_pando_projection"), 1L)
})
