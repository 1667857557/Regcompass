test_that("canonical multitask defaults are method-grounded", {
  defaults <- .rc_multitask_grn_defaults()
  expect_identical(defaults$alpha, 0.5)
  expect_identical(defaults$global_penalty_factor, 1)
  expect_identical(defaults$deviation_penalty_factor, 1)
  expect_identical(defaults$lambda_rule, "lambda.1se")
  expect_identical(defaults$nfolds, 5L)
  expect_identical(defaults$n_bootstrap, 100L)
  expect_identical(defaults$min_selection_frequency, 0.7)
  expect_identical(defaults$min_sign_stability, 0.8)
  expect_identical(defaults$min_detected_cells_per_condition, 10L)
  expect_identical(defaults$min_detection_fraction_per_condition, 0.01)
  expect_identical(defaults$max_edges_per_target, Inf)
  expect_identical(
    formals(.rc_run_celltype_multitask_grns_core)$min_cells,
    100L
  )
})

test_that("canonical candidate universes reject arbitrary top-K truncation", {
  expect_error(
    .rc_validate_multitask_grn_args(list(max_edges_per_target = 100L)),
    "finite top-K would be arbitrary",
    fixed = TRUE
  )
  expect_error(
    .rc_validate_canonical_pando_design_args(list(
      max_edges_per_target = 100L
    )),
    "arbitrary deterministic top-K",
    fixed = TRUE
  )
})

test_that("canonical Pando structural defaults are explicit", {
  defaults <- .rc_validate_canonical_pando_design_args()
  expect_identical(defaults$peak_to_gene_method, "Signac")
  expect_identical(defaults$upstream, 100000)
  expect_identical(defaults$downstream, 0)
  expect_identical(defaults$extend, 1000000)
  expect_false(defaults$only_tss)
  expect_identical(defaults$min_tf_detection, 0)
  expect_identical(defaults$min_peak_detection, 0)
  expect_identical(defaults$min_target_detection, 0)
  expect_identical(defaults$max_edges_per_target, Inf)
})

test_that("candidate observability requires a non-zero TF x peak predictor", {
  cells <- paste0("c", seq_len(40))
  rna <- matrix(
    0,
    nrow = 6,
    ncol = 40,
    dimnames = list(c("TF1", "TF2", "TF3", "G1", "G2", "G3"), cells)
  )
  atac <- matrix(
    0,
    nrow = 3,
    ncol = 40,
    dimnames = list(c("P1", "P2", "P3"), cells)
  )
  meta <- data.frame(
    condition = rep(c("A", "B"), each = 20),
    row.names = cells,
    stringsAsFactors = FALSE
  )

  # E1 has ten same-cell TF/peak observations and ten target observations in A.
  rna["TF1", 1:10] <- 1
  atac["P1", 1:10] <- 1
  rna["G1", 1:10] <- 1

  # E2 is a one-cell event and must be removed.
  rna["TF2", 1] <- 1
  atac["P2", 1] <- 1
  rna["G2", 1] <- 1

  # E3 has ten TF-positive and ten peak-positive cells, but no same-cell overlap.
  rna["TF3", 1:10] <- 1
  atac["P3", 11:20] <- 1
  rna["G3", 1:20] <- 1

  candidates <- data.frame(
    edge_id = c("E1", "E2", "E3"),
    tf_feature_id = c("TF1", "TF2", "TF3"),
    atac_feature_id = c("P1", "P2", "P3"),
    target_feature_id = c("G1", "G2", "G3"),
    target = c("G1", "G2", "G3"),
    stringsAsFactors = FALSE
  )
  args <- .rc_validate_multitask_grn_args()
  filtered <- .rc_filter_shared_candidate_observability(
    candidates = candidates,
    rna = Matrix::Matrix(rna, sparse = TRUE),
    atac = Matrix::Matrix(atac, sparse = TRUE),
    meta = meta,
    condition_col = "condition",
    args = args,
    design_fingerprint = "md5:test"
  )

  expect_identical(filtered$candidates$edge_id, "E1")
  expect_identical(filtered$candidates$n_observable_conditions, 1)
  expect_identical(filtered$candidates$observable_conditions, "A")
  expect_identical(filtered$policy$required_cells_by_condition[["A"]], 10L)
  expect_match(filtered$model_edge_universe_id, "^md5:")
})

test_that("active-edge CV gate must outperform the centred null", {
  expect_identical(
    .rc_cv_predictive_gate(c(-0.1, 0, 0.01, NA), 0),
    c(FALSE, FALSE, TRUE, FALSE)
  )
  expect_identical(
    .rc_cv_predictive_gate(c(0.01, 0.05, 0.1), 0.05),
    c(FALSE, TRUE, TRUE)
  )
})

test_that("bootstrap precision diagnostics use binomial uncertainty", {
  interval <- .rc_binomial_wilson_interval(70, 100)
  expect_equal(interval$lower, 0.604, tolerance = 0.002)
  expect_equal(interval$upper, 0.781, tolerance = 0.002)
  expect_equal((1 + 0.8) / 2, 0.9)
})
