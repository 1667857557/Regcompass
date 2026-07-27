test_that("compact active-edge output retains bootstrap sampling provenance", {
  edges <- data.frame(
    edge_id = "E1",
    condition = "A",
    cell_type = "T",
    tf = "TF1",
    atac_feature_id = "peak1",
    region = "chr1:1-100",
    target = "G1",
    effective_estimate = 0.5,
    selection_frequency = 0.9,
    sign_stability = 1,
    stable_estimate = 0.45,
    cv_rsq = 0.2,
    bootstrap_method = "condition_stratified_sample_cluster_nonparametric",
    bootstrap_resampling_unit = "sample",
    bootstrap_sample_col = "donor",
    n_bootstrap_samples_total = 6L,
    min_bootstrap_samples_per_condition = 3L,
    bootstrap_fallback_reason = NA_character_,
    bootstrap_success_fraction = 1,
    active_edge = TRUE,
    evidence_type = "direct_theta_bootstrap_stability_selected",
    stringsAsFactors = FALSE
  )
  compact <- .rc_compact_active_edges(
    list(tf_peak_gene_significant = edges),
    condition_col = "condition",
    celltype_col = "cell_type"
  )
  expect_true(all(c(
    "bootstrap_method", "bootstrap_resampling_unit", "bootstrap_sample_col",
    "n_bootstrap_samples_total", "min_bootstrap_samples_per_condition",
    "bootstrap_fallback_reason"
  ) %in% colnames(compact)))
  expect_identical(compact$bootstrap_resampling_unit, "sample")
  expect_identical(compact$bootstrap_sample_col, "donor")
})

test_that("sample-aware Stage 1 status is written before downstream core mapping", {
  text <- paste(deparse(body(.rc_run_celltype_multitask_grns)), collapse = "\n")
  expect_match(text, "pando_celltype_status.tsv.gz", fixed = TRUE)
  expect_match(text, "pando_group_status.tsv.gz", fixed = TRUE)
  expect_match(text, "bootstrap_stability_diagnostics.tsv.gz", fixed = TRUE)
  expect_match(text, "bootstrap_policy", fixed = TRUE)
})

test_that("bootstrap policy reaches compact result parameters and provenance", {
  text <- paste(deparse(body(rc_regcompass_step_results)), collapse = "\n")
  policy <- regexpr(
    "bootstrap_policy <- grn$grn_result$bootstrap_policy",
    text,
    fixed = TRUE
  )[[1L]]
  params <- regexpr(
    "bootstrap_resampling_unit",
    text,
    fixed = TRUE
  )[[1L]]
  provenance <- regexpr(
    "stage_provenance",
    text,
    fixed = TRUE
  )[[1L]]
  expect_true(all(c(policy, params, provenance) > 0L))
  expect_true(policy < params)
  expect_true(params < provenance)
  expect_match(text, "bootstrap_fallback_reason", fixed = TRUE)
})

test_that("public workflow forwards sample_col only to Stage 1 and strips timing", {
  text <- paste(deparse(body(rc_run_regcompass)), collapse = "\n")
  expect_match(text, "pando_args$sample_col <- sample_col", fixed = TRUE)
  expect_match(text, "answer$timing <- NULL", fixed = TRUE)
  expect_match(text, "00_execution_timing.tsv", fixed = TRUE)
  expect_false("sample_col" %in% names(formals(rc_regcompass_step_metacells)))
})
