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

  condition <- rep(c("A", "B"), each = 2L)
  estimate <- c(1, -0.2, -0.7, 0.1)
  shared <- c(0.15, -0.05, 0.15, -0.05)
  screen_pval <- c(0.001, 0.02, 0.01, 0.4)
  screen_padj <- c(0.002, 0.04, 0.02, 0.8)
  screen_significant <- screen_padj < padj_threshold
  final_significant <- screen_significant
  coefficient <- data.frame(
    edge_id = rep(edge$edge_id, 2L),
    target = "G",
    tf = rep(edge$tf, 2L),
    region = rep(edge$region, 2L),
    condition = condition,
    estimate = estimate,
    shared_estimate = shared,
    condition_deviation = estimate - shared,
    std_err = NA_real_,
    statistic = NA_real_,
    pval = NA_real_,
    padj = NA_real_,
    screen_estimate = c(0.9, -0.18, -0.65, 0.08),
    screen_std_err = c(0.1, 0.08, 0.14, 0.1),
    screen_statistic = c(9, -2.25, -4.64, 0.8),
    screen_pval = screen_pval,
    screen_padj = screen_padj,
    screen_estimable = TRUE,
    screen_significant = screen_significant,
    significant = final_significant,
    penalty_effect = ifelse(final_significant, estimate, 0),
    estimable = TRUE,
    zero_variance = FALSE,
    aliased = FALSE,
    direction = ifelse(estimate > 0, "positive", "negative"),
    stringsAsFactors = FALSE
  )

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
    contrast_se = NA_real_,
    contrast_statistic = NA_real_,
    contrast_pval = NA_real_,
    contrast_padj = NA_real_,
    contrast_estimable = TRUE,
    contrast_significant = FALSE,
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
    fit_engine = "screened_dictionary_multitask_ridge_effect_refit",
    coefficient_scale = "raw_tf_atac_interaction_units",
    internal_predictor_scale = "equal_condition_within_condition_rms",
    inference_scope =
      "post_screen_joint_ridge_effect_estimation_only_no_final_inference",
    inference_performed = FALSE,
    screening_inference_scope =
      "approximate_ridge_wald_screen_conditional_on_candidate_dictionary_cv_lambda_and_fusion",
    screening_adjust_method = "BH",
    screening_padj_threshold = padj_threshold,
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
      inference_performed = FALSE,
      stringsAsFactors = FALSE
    ),
    network_names = c(A = "net_A", B = "net_B"),
    padj_threshold = padj_threshold,
    adjust_method = "BH",
    scale = FALSE,
    interaction = ":",
    projection_effect_column = "penalty_effect",
    projection_policy = "screen_bh_supported_refit_ridge_effects",
    fit_dictionary_policy =
      "preliminary_joint_ridge_bh_supported_union_then_effect_only_joint_refit",
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

test_that("screened common dictionary is complete in every condition", {
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
    "preliminary BH support"
  )
})

test_that("final refit contains no second-round coefficient inference", {
  fit <- .strict_fit_fixture()
  expect_invisible(RegCompassR:::.rc_require_pando_condition_grn_fit(fit))

  wrong <- fit
  wrong$coefficients$pval[[1L]] <- 0.01
  expect_error(
    RegCompassR:::.rc_require_pando_condition_grn_fit(wrong),
    "must not contain second-round"
  )
})

test_that("final penalty effect uses screen support and final refit beta", {
  fit <- .strict_fit_fixture()
  expect_invisible(RegCompassR:::.rc_require_pando_condition_grn_fit(fit))

  nonsupported <- which(!fit$coefficients$screen_significant)[[1L]]
  expect_equal(fit$coefficients$penalty_effect[[nonsupported]], 0)

  wrong_effect <- fit
  wrong_effect$coefficients$penalty_effect[[nonsupported]] <-
    wrong_effect$coefficients$estimate[[nonsupported]]
  expect_error(
    RegCompassR:::.rc_require_pando_condition_grn_fit(wrong_effect),
    "supported final-refit ridge"
  )

  wrong_flag <- fit
  wrong_flag$coefficients$screen_significant[[1L]] <- FALSE
  expect_error(
    RegCompassR:::.rc_require_pando_condition_grn_fit(wrong_flag),
    "screen_significant"
  )
})

test_that("screening threshold is configurable above 0.1 with no artificial cap", {
  fit <- .strict_fit_fixture(padj_threshold = 0.2)
  expect_invisible(RegCompassR:::.rc_require_pando_condition_grn_fit(fit))

  invalid <- fit
  invalid$padj_threshold <- 1
  invalid$screening_padj_threshold <- 1
  invalid$dictionary_screening_threshold <- 1
  expect_error(
    RegCompassR:::.rc_require_pando_condition_grn_fit(invalid),
    "incomplete"
  )
})

test_that("final condition contrasts are quantitative without second inference", {
  fit <- .strict_fit_fixture()
  expect_invisible(RegCompassR:::.rc_require_pando_condition_grn_fit(fit))
  expect_equal(
    fit$contrasts$contrast_estimate,
    fit$contrasts$estimate_a - fit$contrasts$estimate_b
  )
  expect_true(all(is.na(fit$contrasts$contrast_pval)))
  expect_true(all(is.na(fit$contrasts$contrast_padj)))

  wrong <- fit
  wrong$contrasts$contrast_padj[[1L]] <- 0.01
  expect_error(
    RegCompassR:::.rc_require_pando_condition_grn_fit(wrong),
    "without second-round inference"
  )
})

test_that("Layer 1 gates final refit effects with preliminary screen support", {
  fit <- .strict_fit_fixture()
  gated <- RegCompassR:::.rc_apply_condition_penalty_gate(fit)
  threshold <- fit$padj_threshold
  expected <- fit$coefficients$screen_significant &
    fit$coefficients$estimable &
    is.finite(fit$coefficients$screen_padj) &
    fit$coefficients$screen_padj < threshold

  expect_identical(gated$coefficients$significant, expected)
  expect_equal(
    gated$coefficients$penalty_effect,
    ifelse(expected, gated$coefficients$estimate, 0)
  )
  expect_true(all(gated$coefficients$fit_status == "ok"))
  expect_identical(
    gated$regcompass_penalty_filter,
    paste0(
      "screen_significant & screen_padj < 0.05 & final estimable & finite ",
      "estimate & fit_status == 'ok'"
    )
  )
  expect_identical(
    gated$regcompass_significance_role,
    "stage1_dictionary_screen_support_gate_on_final_refit_effect"
  )
  expect_invisible(
    RegCompassR:::.rc_require_layer1_condition_grn_fit(gated)
  )
})

test_that("non-estimable final condition edges contribute zero", {
  fit <- .strict_fit_fixture()
  fit$coefficients$estimable[[2L]] <- FALSE
  fit$coefficients$penalty_effect[[2L]] <- 0
  fit$coefficients$significant[[2L]] <- FALSE

  gated <- RegCompassR:::.rc_apply_condition_penalty_gate(fit)
  expect_equal(gated$coefficients$penalty_effect[[2L]], 0)
  expect_false(gated$coefficients$significant[[2L]])
  expect_invisible(
    RegCompassR:::.rc_require_layer1_condition_grn_fit(gated)
  )
})

test_that("target reliability requires at least one supported active edge", {
  fit <- .strict_fit_fixture()
  reliability <- RegCompassR:::.rc_condition_target_reliability(fit)
  expect_true(all(reliability$n_significant_edges > 0L))
  expect_true(all(is.finite(reliability$reliability)))

  none <- fit
  b <- none$coefficients$condition == "B"
  none$coefficients$screen_significant[b] <- FALSE
  none$coefficients$screen_padj[b] <- 1
  none$coefficients$significant[b] <- FALSE
  none$coefficients$penalty_effect[b] <- 0
  reliability_none <- RegCompassR:::.rc_condition_target_reliability(none)
  expect_true(is.na(reliability_none$reliability[
    reliability_none$condition == "B"
  ]))
})

test_that("Layer 1 keeps product-of-means projection with screen support gate", {
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
