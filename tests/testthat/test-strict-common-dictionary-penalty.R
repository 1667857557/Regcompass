.strict_fit_fixture <- function(padj_threshold = 0.05, rsq = c(0.8, 0.7)) {
  edge <- data.frame(
    edge_id = c("G||TF1||P1", "G||TF2||P2", "G||TF3||P3"),
    target = "G", tf = c("TF1", "TF2", "TF3"),
    region = c("P1", "P2", "P3"),
    atac_feature_id = c("A1", "A2", "A3"),
    source_global = c(TRUE, TRUE, FALSE),
    source_conditions = c("A;B", "", "A"),
    n_sources = c(3L, 1L, 1L), candidate_index = 1:3,
    stringsAsFactors = FALSE
  )
  support <- data.frame(
    edge_id = c(
      "G||TF1||P1", "G||TF2||P2",
      "G||TF1||P1", "G||TF3||P3", "G||TF1||P1"
    ),
    source_type = c("global", "global", "condition", "condition", "condition"),
    condition = c(NA, NA, "A", "A", "B"),
    peak_target_cor = c(0.4, 0.3, 0.35, 0.25, 0.32),
    tf_target_cor = c(0.5, 0.4, 0.45, 0.35, 0.42),
    stringsAsFactors = FALSE
  )
  support_summary <- data.frame(
    edge_id = edge$edge_id,
    source_global = edge$source_global,
    n_support_conditions = c(2L, 0L, 1L),
    support_conditions = edge$source_conditions,
    n_sources = edge$n_sources,
    stringsAsFactors = FALSE
  )

  condition <- rep(c("A", "B"), each = 3L)
  estimate <- c(1, -0.2, 0.4, 0.8, 0.3, 0.4)
  shared <- rep(colMeans(rbind(estimate[1:3], estimate[4:6])), 2L)
  deviation <- estimate - shared
  pval <- c(0.001, 0.2, 0.5, 0.2, 0.001, 0.5)
  padj <- numeric(6)
  padj[1:3] <- stats::p.adjust(pval[1:3], method = "BH")
  padj[4:6] <- stats::p.adjust(pval[4:6], method = "BH")
  condition_significant <- padj < padj_threshold
  contrast_identifiable <- rep(c(TRUE, TRUE, FALSE), 2L)
  shared_by_boundary <- !contrast_identifiable
  fused_by_penalty <- rep(c(FALSE, FALSE, FALSE), 2L)
  union_edge <- rep(c(TRUE, TRUE, FALSE), 2L)
  supporting <- rep(c("A", "B", ""), 2L)
  n_supporting <- rep(c(1L, 1L, 0L), 2L)

  coefficient <- data.frame(
    edge_id = rep(edge$edge_id, 2L), target = "G",
    tf = rep(edge$tf, 2L), region = rep(edge$region, 2L),
    condition = condition, z = 0.25,
    estimate = estimate, penalty_effect = estimate,
    estimate_standardized = estimate,
    beta_shared = shared, shared_estimate = shared,
    condition_deviation = deviation, delta_beta = deviation,
    contrast_identifiable = contrast_identifiable,
    shared_by_boundary = shared_by_boundary,
    boundary_condition = shared_by_boundary,
    fused_by_penalty = fused_by_penalty,
    fusion_component_id = rep(c("component1", "component1", "component1"), 2L),
    shared_edge = rep(c(FALSE, FALSE, TRUE), 2L),
    raw_information_condition = c(30, 20, 0, 15, 10, 0),
    profile_information = rep(c(10, 6, 0), 2L),
    profile_information_delta = rep(c(10, 6, 0), 2L),
    inference_schema = "scheme_e_fusion_component_joint_refit_v1",
    inference_component_id = rep(c("edge1_component1", "edge2_component1", "edge3_component1"), 2L),
    inference_hypothesis_id = rep(c("H0_e1", "H0_e2", "H0_e3"), 2L),
    inference_estimate = estimate,
    inference_se = rep(0.1, 6L),
    inference_statistic = estimate / 0.1,
    inference_estimable = TRUE,
    std_err = rep(0.1, 6L), statistic = estimate / 0.1,
    pval = pval, padj = padj,
    bh_scope = "condition_target_BH", bh_family_size = 3L,
    condition_significant = condition_significant,
    statistically_supported = condition_significant,
    significant = condition_significant,
    pando_estimation_active = TRUE, active = TRUE,
    edge_union_supported = union_edge,
    supporting_conditions = supporting,
    n_supporting_conditions = n_supporting,
    all_conditions_fit_valid = TRUE,
    active_in_regcompass = union_edge,
    fit_status = "ok",
    estimable = TRUE, zero_variance = FALSE,
    condition_informative = TRUE,
    global_support = rep(edge$source_global, 2L),
    local_support = c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE),
    penalty_family = "information_scaled_sparse_deviation",
    penalty_value = 0.25, solver_status = "ok",
    objective = 1, kkt_residual = 1e-10, iterations = 24L,
    stringsAsFactors = FALSE
  )

  contrast_estimate <- estimate[4:6] - estimate[1:3]
  contrasts <- data.frame(
    edge_id = edge$edge_id, target = "G", tf = edge$tf,
    region = edge$region, atac_feature_id = edge$atac_feature_id,
    condition_a = "A", condition_b = "B", contrast = "B-A",
    estimate_a = estimate[1:3], estimate_b = estimate[4:6],
    contrast_estimate = contrast_estimate,
    inference_contrast_estimate = contrast_estimate,
    contrast_se = c(0.2, 0.2, NA),
    contrast_statistic = c(-1, 2.5, NA),
    contrast_pval = c(0.3, 0.01, NA),
    contrast_padj = c(0.3, 0.02, NA),
    contrast_estimable = c(TRUE, TRUE, FALSE),
    contrast_identifiable = c(TRUE, TRUE, FALSE),
    contrast_status = c("ok", "ok", "fused_by_E"),
    contrast_significant = c(FALSE, TRUE, FALSE),
    shared_by_boundary = c(FALSE, FALSE, TRUE),
    fused_by_penalty = FALSE,
    shared_edge = c(FALSE, FALSE, TRUE),
    profile_information_delta = c(10, 6, 0),
    penalty_family = "information_scaled_sparse_deviation",
    penalty_value = 0.25, solver_status = "ok",
    kkt_residual = 1e-10, iterations = 24L,
    stringsAsFactors = FALSE
  )

  structure(list(
    schema_version = "pando_condition_grn_common_dictionary_v1",
    model_schema = "pando_condition_grn_Estar_jointse_v1",
    fit_engine = "condition_union_Estar_z025_jointse",
    inference_schema = "scheme_e_fusion_component_joint_refit_v1",
    coefficient_scale = "raw_tf_atac_interaction_units",
    internal_predictor_scale = "equal_condition_within_condition_rms",
    inference_scope = "E_star_z025_primary_fusion_component_joint_refit_condition_target_BH",
    cell_type = "T_cell", condition_levels = c("A", "B"),
    reference_condition = "A", condition_col = "condition",
    cell_type_col = "cell_type",
    condition_cell_ids = list(
      A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3")
    ),
    edge_dictionary = edge,
    dictionary_support_table = support,
    dictionary_support_summary = support_summary,
    candidate_tf_cor = 0.1, candidate_peak_cor = 0.05,
    coefficients = coefficient, contrasts = contrasts,
    fit = data.frame(
      target = rep("G", 2), condition = c("A", "B"),
      rsq = rsq, rsq_definition = "scheme_e_z025_full_data_R2_diagnostic",
      fit_status = "ok", sigma2_common = 0.2,
      inference_sigma2 = 0.2, inference_residual_df = 4L,
      inference_rank = 4L, reference_condition = "A",
      deviation_z = 0.25,
      penalty_family = "information_scaled_sparse_deviation",
      solver_status = "ok", kkt_residual = 1e-10, iterations = 24L,
      predictor_scale_reference = "equal_condition_within_condition_rms",
      inference_schema = "scheme_e_fusion_component_joint_refit_v1",
      orthogonality_error = 1e-12, dr_error = 1e-12,
      stringsAsFactors = FALSE
    ),
    network_names = c(A = "net_A", B = "net_B"),
    padj_threshold = padj_threshold, adjust_method = "BH",
    scale = FALSE, interaction = ":",
    projection_effect_column = "penalty_effect",
    projection_policy = "any_condition_padj_exact_edge_union",
    fit_dictionary_policy =
      "global_and_condition_union_pando_correlation_supported_frozen_dictionary",
    candidate_edge_count = 3L, fit_dictionary_edge_count = 3L,
    rna_layer = "data", peak_layer = "data",
    peak_value_type = "normalized",
    preprocessing_fingerprint = "fixture-preprocessing",
    target_genes = "G",
    deviation_penalty = list(
      family = "information_scaled_sparse_deviation", z = 0.25
    ),
    target_solver = list(list(
      status = "ok", kkt_residual = 1e-10, iterations = 24L,
      penalty_family = "information_scaled_sparse_deviation",
      penalty_value = 0.25
    )),
    target_scaling = list(list(
      reference = "equal_condition_within_condition_rms"
    )),
    target_contrast_tree = list(data.frame()),
    rsq_definition = "scheme_e_z025_full_data_R2_diagnostic"
  ), class = c("ConditionGRNFit", "list"))
}

test_that("condition fit contract requires fixed E-star/JSE z=0.25", {
  fit <- .strict_fit_fixture()
  expect_false("lambda" %in% colnames(fit$fit))
  expect_false("rsq_oof" %in% colnames(fit$fit))
  expect_equal(fit$deviation_penalty$z, 0.25)
  expect_invisible(RegCompassR:::.rc_require_pando_condition_grn_fit(fit))
  wrong <- fit
  wrong$deviation_penalty$z <- 0.5
  expect_error(RegCompassR:::.rc_require_pando_condition_grn_fit(wrong),
               "z=0.25")
})

test_that("deduplicated exact dictionary is complete in every condition", {
  fit <- .strict_fit_fixture()
  expect_invisible(RegCompassR:::.rc_require_pando_condition_grn_fit(fit))
  incomplete <- fit
  incomplete$coefficients <- incomplete$coefficients[-1L, , drop = FALSE]
  expect_error(RegCompassR:::.rc_require_pando_condition_grn_fit(incomplete),
               "every condition")
})

test_that("condition-target BH controls union admission without zeroing beta_E", {
  fit <- .strict_fit_fixture()
  expect_invisible(RegCompassR:::.rc_require_pando_condition_grn_fit(fit))
  nonsig_b_e1 <- which(
    fit$coefficients$condition == "B" &
      fit$coefficients$edge_id == "G||TF1||P1"
  )
  expect_false(fit$coefficients$condition_significant[[nonsig_b_e1]])
  expect_true(fit$coefficients$active[[nonsig_b_e1]])
  expect_true(fit$coefficients$active_in_regcompass[[nonsig_b_e1]])
  expect_equal(fit$coefficients$penalty_effect[[nonsig_b_e1]],
               fit$coefficients$estimate[[nonsig_b_e1]])
  e3 <- fit$coefficients$edge_id == "G||TF3||P3"
  expect_false(any(fit$coefficients$active_in_regcompass[e3]))
})

test_that("RegCompass handoff keeps the same admitted exact edges in every condition", {
  fit <- .strict_fit_fixture(rsq = c(0.8, 0.01))
  gated <- RegCompassR:::.rc_apply_condition_penalty_gate(fit)
  expect_identical(gated$coefficients$active, rep(TRUE, 6L))
  expect_equal(gated$coefficients$penalty_effect, fit$coefficients$estimate)
  rows_b <- gated$coefficients$condition == "B"
  expect_false(any(gated$coefficients$target_model_supported[rows_b]))
  expect_identical(
    gated$coefficients$penalty_eligible,
    gated$coefficients$active_in_regcompass
  )
  admitted_a <- gated$coefficients$edge_id[
    gated$coefficients$condition == "A" & gated$coefficients$penalty_eligible
  ]
  admitted_b <- gated$coefficients$edge_id[
    gated$coefficients$condition == "B" & gated$coefficients$penalty_eligible
  ]
  expect_setequal(admitted_a, admitted_b)
  expect_identical(
    gated$regcompass_target_rsq_definition,
    "scheme_e_z025_full_data_R2_diagnostic"
  )
  expect_invisible(RegCompassR:::.rc_require_layer1_condition_grn_fit(gated))
})

test_that("target reliability follows admitted union edges, not target R2", {
  fit <- .strict_fit_fixture(rsq = c(0.8, 0.01))
  gated <- RegCompassR:::.rc_apply_condition_penalty_gate(fit)
  reliability <- RegCompassR:::.rc_condition_target_reliability(gated)
  expect_true(all(reliability$n_projection_edges > 0L))
  expect_equal(reliability$reliability, c(1, 1))
  expect_equal(reliability$target_rsq_supported_diagnostic, c(TRUE, FALSE))
})

test_that("pairwise contrasts close on the same joint production coefficients", {
  fit <- .strict_fit_fixture()
  expect_equal(
    fit$contrasts$contrast_estimate,
    fit$contrasts$estimate_b - fit$contrasts$estimate_a
  )
  expect_false(fit$contrasts$contrast_identifiable[[3L]])
  expect_true(fit$contrasts$shared_by_boundary[[3L]])
  expect_invisible(RegCompassR:::.rc_require_pando_condition_grn_fit(fit))
})

test_that("Layer 1 uses admitted union beta_E with unchanged exposure mixture", {
  selector <- paste(
    deparse(body(RegCompassR:::.rc_condition_pando_object_for_fit)),
    collapse = "\n"
  )
  projection <- paste(
    deparse(body(RegCompassR:::.rc_condition_pando_projection)),
    collapse = "\n"
  )
  helper <- paste(
    deparse(body(RegCompassR:::.rc_pando_projection_from_group_means)),
    collapse = "\n"
  )
  expect_match(selector, ".rc_require_layer1_condition_grn_fit", fixed = TRUE)
  expect_match(projection, ".rc_condition_penalty_gate", fixed = TRUE)
  expect_match(
    projection,
    "beta_times_0.75_paired_mean_product_plus_0.25_product_of_means",
    fixed = TRUE
  )
  expect_match(helper, "paired_product <- tf_block * peak_block", fixed = TRUE)
  expect_match(helper, "product_of_means_weight * tf_mean * peak_mean", fixed = TRUE)
  expect_identical(
    RegCompassR:::.RC_PANDO_PROJECTION_PRODUCT_OF_MEANS_WEIGHT, 0.25
  )
})

test_that("condition contract is implemented once in functional source files", {
  root <- testthat::test_path("..", "..")
  r_files <- list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)
  source_text <- unlist(lapply(r_files, readLines, warn = FALSE), use.names = FALSE)
  count_definition <- function(name) {
    pattern <- paste0(
      "^", gsub("\\.", "\\\\.", name),
      "[[:space:]]*<-[[:space:]]*function"
    )
    sum(grepl(pattern, source_text))
  }
  expect_identical(count_definition(".rc_require_pando_condition_grn_fit"), 1L)
  expect_identical(count_definition(".rc_require_layer1_condition_grn_fit"), 1L)
  expect_identical(count_definition(".rc_extract_condition_grn_contract"), 1L)
  expect_identical(count_definition(".rc_condition_pando_projection"), 1L)
  expect_identical(count_definition(".rc_cell_first_projection_layer1"), 1L)
  expect_identical(count_definition(".rc_cell_first_projection_layer1_v6"), 0L)
})
