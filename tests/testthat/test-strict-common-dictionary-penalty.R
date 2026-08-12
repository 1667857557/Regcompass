.strict_fit_fixture <- function(padj_threshold = 0.05) {
  edge <- data.frame(
    edge_id = c("G||TF1||P1", "G||TF2||P2", "G||TF3||P3"),
    target = "G",
    tf = c("TF1", "TF2", "TF3"),
    region = c("P1", "P2", "P3"),
    atac_feature_id = c("A1", "A2", "A3"),
    source_global = c(TRUE, TRUE, FALSE),
    source_conditions = c("A;B", "", "A"),
    n_sources = c(3L, 1L, 1L),
    max_abs_peak_target_cor = c(0.4, 0.3, 0.25),
    max_abs_tf_target_cor = c(0.5, 0.4, 0.35),
    candidate_index = 1:3,
    stringsAsFactors = FALSE
  )
  attr(edge, "preprocessing_provenance_verified") <- TRUE

  support <- data.frame(
    edge_id = c(
      "G||TF1||P1", "G||TF2||P2",
      "G||TF1||P1", "G||TF3||P3", "G||TF1||P1"
    ),
    target = "G",
    tf = c("TF1", "TF2", "TF1", "TF3", "TF1"),
    region = c("P1", "P2", "P1", "P3", "P1"),
    atac_feature_id = c("A1", "A2", "A1", "A3", "A1"),
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
    max_abs_peak_target_cor = edge$max_abs_peak_target_cor,
    max_abs_tf_target_cor = edge$max_abs_tf_target_cor,
    stringsAsFactors = FALSE
  )

  condition <- rep(c("A", "B"), each = 3L)
  pval <- c(0.001, 0.02, 0.03, 0.001, 0.02, 0.03)
  padj <- unlist(lapply(split(pval, condition), stats::p.adjust,
                        method = "BH"), use.names = FALSE)
  estimate <- c(1, -0.2, 0.4, -0.7, 0.3, 0.6)
  shared <- rep(colMeans(rbind(estimate[1:3], estimate[4:6])), 2L)
  statistically_supported <- padj < padj_threshold
  global_support <- rep(edge$source_global, 2L)
  local_support <- c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE)
  active <- statistically_supported & (global_support | local_support)
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
    statistically_supported = statistically_supported,
    global_support = global_support,
    local_support = local_support,
    active = active,
    significant = active,
    penalty_effect = ifelse(active, estimate, 0),
    estimable = TRUE,
    zero_variance = FALSE,
    aliased = FALSE,
    direction = ifelse(estimate > 0, "positive", "negative"),
    stringsAsFactors = FALSE
  )

  contrast_pval <- c(0.01, 0.2, 0.5)
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
    estimate_a = estimate[1:3],
    estimate_b = estimate[4:6],
    contrast_estimate = estimate[1:3] - estimate[4:6],
    contrast_se = c(0.2, 0.2, 0.2),
    contrast_statistic = c(8.5, -2.5, -1),
    contrast_pval = contrast_pval,
    contrast_padj = contrast_padj,
    contrast_estimable = TRUE,
    contrast_significant = contrast_padj < padj_threshold,
    stringsAsFactors = FALSE
  )

  structure(list(
    schema_version = "pando_condition_grn_common_dictionary_v1",
    model_schema = "pando_condition_grn_multitask_ridge_v3",
    fit_engine = "condition_union_single_no_fusion_common_lambda_ridge",
    coefficient_scale = "raw_tf_atac_interaction_units",
    internal_predictor_scale = "equal_condition_within_condition_rms",
    inference_scope =
      "approximate_ridge_wald_conditional_on_global_or_condition_pando_screened_dictionary_and_cv_lambda",
    cell_type = "T_cell",
    condition_levels = c("A", "B"),
    condition_col = "condition",
    cell_type_col = "cell_type",
    condition_cell_ids = list(
      A = c("a1", "a2", "a3"),
      B = c("b1", "b2", "b3")
    ),
    edge_dictionary = edge,
    dictionary_support_table = support,
    dictionary_support_summary = support_summary,
    candidate_tf_cor = 0.1,
    candidate_peak_cor = 0.05,
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
    projection_policy =
      "active_global_or_local_pando_support_and_condition_bh_ridge_effects",
    fit_dictionary_policy =
      "global_and_condition_union_pando_correlation_supported_frozen_dictionary",
    candidate_edge_count = 3L,
    fit_dictionary_edge_count = 3L,
    rna_layer = "data",
    peak_layer = "data",
    peak_value_type = "normalized",
    preprocessing_fingerprint = "fixture-preprocessing",
    target_genes = "G"
  ), class = c("ConditionGRNFit", "list"))
}

test_that("deduplicated common dictionary is complete in every condition", {
  fit <- .strict_fit_fixture()
  expect_invisible(RegCompassR:::.rc_require_pando_condition_grn_fit(fit))
  incomplete <- fit
  incomplete$coefficients <- incomplete$coefficients[-1L, , drop = FALSE]
  expect_error(
    RegCompassR:::.rc_require_pando_condition_grn_fit(incomplete),
    "every common-dictionary edge"
  )
  duplicated <- fit
  duplicated$edge_dictionary$edge_id[[2L]] <- duplicated$edge_dictionary$edge_id[[1L]]
  expect_error(
    RegCompassR:::.rc_require_pando_condition_grn_fit(duplicated),
    "common dictionary"
  )
})

test_that("global-only support can rescue a small-condition candidate miss", {
  fit <- .strict_fit_fixture()
  expect_invisible(RegCompassR:::.rc_require_pando_condition_grn_fit(fit))
  row <- which(fit$coefficients$condition == "B" &
               fit$coefficients$edge_id == "G||TF2||P2")
  expect_true(fit$coefficients$statistically_supported[[row]])
  expect_true(fit$coefficients$global_support[[row]])
  expect_false(fit$coefficients$local_support[[row]])
  expect_true(fit$coefficients$active[[row]])
})

test_that("local-only edge is not activated in an unsupported condition", {
  fit <- .strict_fit_fixture()
  row <- which(fit$coefficients$condition == "B" &
               fit$coefficients$edge_id == "G||TF3||P3")
  expect_true(fit$coefficients$statistically_supported[[row]])
  expect_false(fit$coefficients$global_support[[row]])
  expect_false(fit$coefficients$local_support[[row]])
  expect_false(fit$coefficients$active[[row]])
  expect_equal(fit$coefficients$penalty_effect[[row]], 0)
})

test_that("condition activity uses configured BH threshold plus dictionary support", {
  fit <- .strict_fit_fixture()
  expect_invisible(RegCompassR:::.rc_require_pando_condition_grn_fit(fit))

  wrong_effect <- fit
  row <- which(!fit$coefficients$active)[[1L]]
  wrong_effect$coefficients$penalty_effect[[row]] <-
    wrong_effect$coefficients$estimate[[row]]
  expect_error(
    RegCompassR:::.rc_require_pando_condition_grn_fit(wrong_effect),
    "penalty_effect"
  )

  wrong_flag <- fit
  wrong_flag$coefficients$active[[1L]] <- FALSE
  wrong_flag$coefficients$significant[[1L]] <- FALSE
  expect_error(
    RegCompassR:::.rc_require_pando_condition_grn_fit(wrong_flag),
    "activity flags"
  )
})

test_that("padj threshold remains configurable over the valid interval", {
  fit <- .strict_fit_fixture(padj_threshold = 0.2)
  expect_invisible(RegCompassR:::.rc_require_pando_condition_grn_fit(fit))
  invalid <- fit
  invalid$padj_threshold <- 1
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

test_that("Layer 1 preserves Pando activity and adds only fit-status eligibility", {
  fit <- .strict_fit_fixture()
  gated <- RegCompassR:::.rc_apply_condition_penalty_gate(fit)
  expect_identical(gated$coefficients$active, fit$coefficients$active)
  expect_identical(gated$coefficients$significant, fit$coefficients$significant)
  expect_equal(gated$coefficients$penalty_effect, fit$coefficients$penalty_effect)
  expect_identical(
    gated$coefficients$penalty_eligible,
    fit$coefficients$active
  )
  expect_true(all(gated$coefficients$fit_status == "ok"))
  expect_identical(
    gated$regcompass_penalty_filter,
    "Pando active edge & target fit_status == 'ok'"
  )
  expect_identical(
    gated$regcompass_significance_role,
    "consume_pando_active_condition_edge_without_reselection"
  )
  expect_invisible(
    RegCompassR:::.rc_require_layer1_condition_grn_fit(gated)
  )
})

test_that("non-estimable condition edges remain inactive", {
  fit <- .strict_fit_fixture()
  fit$coefficients$estimable[[2L]] <- FALSE
  fit$coefficients$statistically_supported[[2L]] <- FALSE
  fit$coefficients$active[[2L]] <- FALSE
  fit$coefficients$significant[[2L]] <- FALSE
  fit$coefficients$penalty_effect[[2L]] <- 0
  fit$coefficients$padj[[2L]] <- NA_real_
  fit$coefficients$pval[[2L]] <- NA_real_

  gated <- RegCompassR:::.rc_apply_condition_penalty_gate(fit)
  expect_equal(gated$coefficients$penalty_effect[[2L]], 0)
  expect_false(gated$coefficients$active[[2L]])
  expect_false(gated$coefficients$penalty_eligible[[2L]])
  expect_invisible(
    RegCompassR:::.rc_require_layer1_condition_grn_fit(gated)
  )
})

test_that("target reliability requires at least one active edge", {
  fit <- .strict_fit_fixture()
  reliability <- RegCompassR:::.rc_condition_target_reliability(fit)
  expect_true(all(reliability$n_projection_edges > 0L))
  expect_true(all(is.finite(reliability$reliability)))

  none <- fit
  rows <- none$coefficients$condition == "B"
  none$coefficients$padj[rows] <- 1
  none$coefficients$statistically_supported[rows] <- FALSE
  none$coefficients$active[rows] <- FALSE
  none$coefficients$significant[rows] <- FALSE
  none$coefficients$penalty_effect[rows] <- 0
  reliability_none <- RegCompassR:::.rc_condition_target_reliability(none)
  expect_true(is.na(reliability_none$reliability[
    reliability_none$condition == "B"
  ]))
})

test_that("Layer 1 keeps product-of-means projection with Pando active gate", {
  selector <- paste(
    deparse(body(RegCompassR:::.rc_condition_pando_object_for_fit)),
    collapse = "\n"
  )
  projection <- paste(
    deparse(body(RegCompassR:::.rc_condition_pando_projection)),
    collapse = "\n"
  )
  expect_match(selector, ".rc_require_layer1_condition_grn_fit", fixed = TRUE)
  expect_match(projection, ".rc_require_layer1_condition_grn_fit", fixed = TRUE)
  expect_match(projection, ".rc_condition_penalty_gate", fixed = TRUE)
  expect_match(
    projection, "beta_times_group_mean_tf_times_group_mean_atac", fixed = TRUE
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
