.strict_fit_fixture <- function(padj_threshold = 0.05, rsq = c(0.8, 0.7)) {
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
  estimate <- c(1, -0.2, 0.4, -0.7, 0.3, 0.4)
  shared_edge <- colMeans(rbind(estimate[1:3], estimate[4:6]))
  shared <- rep(shared_edge, 2L)
  deviation <- estimate - shared
  pval <- c(0.001, 0.2, NA, 0.001, 0.02, NA)
  padj <- unlist(lapply(split(seq_along(pval), condition), function(index) {
    value <- rep(NA_real_, length(index))
    valid <- is.finite(pval[index])
    value[valid] <- stats::p.adjust(pval[index][valid], method = "BH")
    value
  }), use.names = FALSE)
  statistically_supported <- is.finite(padj) & padj < padj_threshold
  global_support <- rep(edge$source_global, 2L)
  local_support <- c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE)
  identifiable_edge <- c(TRUE, TRUE, FALSE)
  contrast_identifiable <- rep(identifiable_edge, 2L)
  shared_by_boundary <- !contrast_identifiable
  fused_by_penalty <- rep(c(FALSE, TRUE, FALSE), 2L)
  coefficient <- data.frame(
    edge_id = rep(edge$edge_id, 2L),
    target = "G",
    tf = rep(edge$tf, 2L),
    region = rep(edge$region, 2L),
    condition = condition,
    estimate = estimate,
    estimate_standardized = estimate,
    shared_estimate = shared,
    beta_shared = shared,
    condition_deviation = deviation,
    delta_beta = deviation,
    std_err = c(0.1, 0.1, NA, 0.1, 0.1, NA),
    statistic = c(10, -2, NA, -7, 3, NA),
    pval = pval,
    padj = padj,
    statistically_supported = statistically_supported,
    global_support = global_support,
    local_support = local_support,
    active = TRUE,
    significant = statistically_supported,
    penalty_effect = estimate,
    estimable = TRUE,
    zero_variance = FALSE,
    condition_informative = TRUE,
    contrast_identifiable = contrast_identifiable,
    shared_by_boundary = shared_by_boundary,
    fused_by_penalty = fused_by_penalty,
    raw_information_condition = c(30, 20, 0, 15, 10, 0),
    profile_information_delta = rep(c(10, 6, 0), 2L),
    profile_information_definition = "pairwise_delta_profile_information",
    penalty_family = "exact_edge_sparse_deviation",
    penalty_value = 0.25,
    solver_status = "ok",
    kkt_residual = 1e-10,
    iterations = 24L,
    aliased = !contrast_identifiable,
    direction = ifelse(estimate > 0, "positive",
                       ifelse(estimate < 0, "negative", "zero")),
    stringsAsFactors = FALSE
  )

  contrast_pval <- c(0.01, 0.2, NA)
  valid_contrast <- is.finite(contrast_pval) & identifiable_edge
  contrast_padj <- rep(NA_real_, 3L)
  contrast_padj[valid_contrast] <- stats::p.adjust(
    contrast_pval[valid_contrast], method = "BH"
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
    estimate_a = estimate[1:3],
    estimate_b = estimate[4:6],
    contrast_estimate = estimate[1:3] - estimate[4:6],
    contrast_se = c(0.2, 0.2, NA),
    contrast_statistic = c(8.5, -2.5, NA),
    contrast_pval = contrast_pval,
    contrast_padj = contrast_padj,
    contrast_estimable = identifiable_edge,
    contrast_identifiable = identifiable_edge,
    shared_by_boundary = !identifiable_edge,
    fused_by_penalty = c(FALSE, TRUE, FALSE),
    profile_information_delta = c(10, 6, 0),
    penalty_family = "exact_edge_sparse_deviation",
    penalty_value = 0.25,
    solver_status = "ok",
    kkt_residual = 1e-10,
    iterations = 24L,
    contrast_significant = identifiable_edge &
      is.finite(contrast_padj) & contrast_padj < padj_threshold,
    stringsAsFactors = FALSE
  )

  structure(list(
    schema_version = "pando_condition_grn_common_dictionary_v1",
    model_schema = "pando_condition_grn_sparse_deviation_v4",
    fit_engine = "condition_union_scheme_e_exact_edge_z025",
    coefficient_scale = "raw_tf_atac_interaction_units",
    internal_predictor_scale = "equal_condition_within_condition_rms",
    inference_scope = "scheme_e_z025_primary;BH_and_R2_are_diagnostics_only",
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
      rsq = rsq, fit_status = "ok", sigma2_common = 0.2,
      deviation_z = 0.25,
      penalty_family = "exact_edge_sparse_deviation",
      solver_status = "ok", kkt_residual = 1e-10, iterations = 24L,
      predictor_scale_reference = "equal_condition_within_condition_rms",
      profile_information_definition = "pairwise_delta_profile_information",
      stringsAsFactors = FALSE
    ),
    network_names = c(A = "net_A", B = "net_B"),
    padj_threshold = padj_threshold,
    adjust_method = "BH",
    scale = FALSE,
    interaction = ":",
    projection_effect_column = "penalty_effect",
    projection_policy = "continuous_common_dictionary_scheme_e_effects",
    fit_dictionary_policy =
      "global_and_condition_union_pando_correlation_supported_frozen_dictionary",
    candidate_edge_count = 3L,
    fit_dictionary_edge_count = 3L,
    rna_layer = "data",
    peak_layer = "data",
    peak_value_type = "normalized",
    preprocessing_fingerprint = "fixture-preprocessing",
    target_genes = "G",
    deviation_penalty = list(
      family = "exact_edge_sparse_deviation", z = 0.25
    ),
    target_solver = list(list(
      status = "ok", kkt_residual = 1e-10, iterations = 24L,
      penalty_family = "exact_edge_sparse_deviation", penalty_value = 0.25
    )),
    target_scaling = list(list(
      reference = "equal_condition_within_condition_rms"
    )),
    rsq_definition = "scheme_e_z025_full_data_R2_diagnostic"
  ), class = c("ConditionGRNFit", "list"))
}

test_that("condition fit contract requires fixed Scheme E z=0.25", {
  fit <- .strict_fit_fixture()
  expect_false("lambda" %in% colnames(fit$fit))
  expect_false("rsq_oof" %in% colnames(fit$fit))
  expect_equal(fit$deviation_penalty$z, 0.25)
  expect_invisible(RegCompassR:::.rc_require_pando_condition_grn_fit(fit))

  wrong <- fit
  wrong$deviation_penalty$z <- 0.5
  expect_error(
    RegCompassR:::.rc_require_pando_condition_grn_fit(wrong),
    "z=0.25"
  )
})

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
  duplicated$edge_dictionary$edge_id[[2L]] <-
    duplicated$edge_dictionary$edge_id[[1L]]
  expect_error(
    RegCompassR:::.rc_require_pando_condition_grn_fit(duplicated),
    "common dictionary"
  )
})

test_that("BH remains diagnostic and cannot change Scheme E topology", {
  fit <- .strict_fit_fixture()
  expect_invisible(RegCompassR:::.rc_require_pando_condition_grn_fit(fit))
  unsupported <- which(!fit$coefficients$statistically_supported)[[1L]]
  expect_true(fit$coefficients$active[[unsupported]])
  expect_false(fit$coefficients$significant[[unsupported]])
  expect_equal(
    fit$coefficients$penalty_effect[[unsupported]],
    fit$coefficients$estimate[[unsupported]]
  )

  wrong <- fit
  wrong$coefficients$active[[unsupported]] <- FALSE
  expect_error(
    RegCompassR:::.rc_require_pando_condition_grn_fit(wrong),
    "finite common-dictionary"
  )
})

test_that("global and local support remain provenance only", {
  fit <- .strict_fit_fixture()
  global_only <- which(
    fit$coefficients$condition == "B" &
      fit$coefficients$edge_id == "G||TF2||P2"
  )
  local_admitted <- which(
    fit$coefficients$condition == "B" &
      fit$coefficients$edge_id == "G||TF3||P3"
  )
  expect_true(fit$coefficients$global_support[[global_only]])
  expect_false(fit$coefficients$local_support[[global_only]])
  expect_true(fit$coefficients$active[[global_only]])
  expect_false(fit$coefficients$global_support[[local_admitted]])
  expect_false(fit$coefficients$local_support[[local_admitted]])
  expect_true(fit$coefficients$active[[local_admitted]])
})

test_that("zero-information contrast is exact-shared and not claim-eligible", {
  fit <- .strict_fit_fixture()
  boundary <- fit$coefficients$shared_by_boundary %in% TRUE
  expect_true(all(!fit$coefficients$contrast_identifiable[boundary]))
  expect_equal(fit$coefficients$condition_deviation[boundary], c(0, 0))
  expect_equal(
    fit$contrasts$contrast_estimate[fit$contrasts$shared_by_boundary], 0
  )
  wrong <- fit
  index <- which(wrong$coefficients$shared_by_boundary)[[1L]]
  wrong$coefficients$condition_deviation[[index]] <- 0.1
  expect_error(
    RegCompassR:::.rc_require_pando_condition_grn_fit(wrong),
    "exact-shared"
  )
})

test_that("RegCompass preserves continuous effects and makes target R2 diagnostic only", {
  fit <- .strict_fit_fixture(rsq = c(0.8, 0.01))
  gated <- RegCompassR:::.rc_apply_condition_penalty_gate(fit)
  expect_identical(gated$coefficients$active, rep(TRUE, 6L))
  expect_equal(gated$coefficients$penalty_effect, fit$coefficients$estimate)
  rows_b <- gated$coefficients$condition == "B"
  expect_false(any(gated$coefficients$target_model_supported[rows_b]))
  expect_true(all(gated$coefficients$penalty_eligible[rows_b]))
  expect_identical(
    gated$regcompass_penalty_filter,
    "finite continuous Scheme-E coefficient on frozen dictionary & fit_status == 'ok'"
  )
  expect_identical(
    gated$regcompass_target_rsq_definition,
    "scheme_e_z025_full_data_R2_diagnostic"
  )
  expect_invisible(RegCompassR:::.rc_require_layer1_condition_grn_fit(gated))
})

test_that("target reliability is availability, not BH or R2 filtering", {
  fit <- .strict_fit_fixture(rsq = c(0.8, 0.01))
  reliability <- RegCompassR:::.rc_condition_target_reliability(fit)
  expect_true(all(reliability$n_projection_edges > 0L))
  expect_equal(reliability$reliability, c(1, 1))
  expect_equal(
    reliability$target_rsq_supported_diagnostic,
    c(TRUE, FALSE)
  )

  no_bh <- fit
  rows <- no_bh$coefficients$condition == "B"
  no_bh$coefficients$padj[rows] <- 1
  no_bh$coefficients$statistically_supported[rows] <- FALSE
  no_bh$coefficients$significant[rows] <- FALSE
  reliability_no_bh <- RegCompassR:::.rc_condition_target_reliability(no_bh)
  expect_equal(
    reliability_no_bh$reliability[reliability_no_bh$condition == "B"], 1
  )
})

test_that("pairwise differential claims require identifiable Scheme E contrasts", {
  fit <- .strict_fit_fixture()
  expect_equal(
    fit$contrasts$contrast_estimate,
    fit$contrasts$estimate_a - fit$contrasts$estimate_b
  )
  expect_false(fit$contrasts$contrast_identifiable[[3L]])
  expect_true(fit$contrasts$shared_by_boundary[[3L]])
  wrong <- fit
  wrong$contrasts$contrast_padj[[1L]] <- 0.99
  expect_error(
    RegCompassR:::.rc_require_pando_condition_grn_fit(wrong),
    "contrast padj"
  )
})

test_that("Layer 1 keeps paired-shrinkage exposure but uses continuous Scheme E beta", {
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
    RegCompassR:::.RC_PANDO_PROJECTION_PRODUCT_OF_MEANS_WEIGHT,
    0.25
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
