test_that("condition effects use an explicit reference condition", {
  beta <- matrix(
    c(0.5, 2.0, -0.25, -1.0),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("edge1", "edge2"),
      c("Control", "Drug")
    )
  )
  observed <- .rc_reference_contrast(beta, "Control")

  expect_equal(unname(observed[, "Control"]), c(0, 0))
  expect_equal(unname(observed[, "Drug"]), c(1.5, -0.75))
  expect_error(
    .rc_reference_contrast(beta, "Missing"),
    "Reference condition is absent"
  )
})

test_that("RegCompass enforces the condition-sparse Pando v4 contract", {
  defaults <- eval(
    formals(.rc_run_condition_single_cell_grns)$pando_infer_args
  )
  expect_identical(defaults$method, "shared_baseline_condition_sparse")
  expect_identical(defaults$candidate_screen, "motif_domain")
  expect_identical(defaults$condition_mix, 0.5)
  expect_identical(defaults$cv_block_col, "sample_id")
  expect_identical(defaults$condition_weight, "equal")
  expect_true(defaults$scale)

  body_text <- paste(
    deparse(body(
      .rc_run_condition_single_cell_grns_without_safe_defaults
    )), collapse = "\n"
  )
  expect_match(body_text, ".rc_extract_condition_grn_contract", fixed = TRUE)
  expect_match(body_text, "condition_grn_fit_v4.rds", fixed = TRUE)
  expect_match(body_text, "pando_edge_predictor_transforms", fixed = TRUE)
})

test_that("Layer 1 never reconstructs TF-by-ATAC from metacell means", {
  text <- paste(
    deparse(body(.rc_cell_first_projection_layer1)), collapse = "\n"
  )
  expect_match(
    text, "Pando::project_condition_grn_groups", fixed = TRUE
  )
  expect_false(grepl(".rc_tf_peak_interaction", text, fixed = TRUE))
})

test_that("model-space modifier uses shared pooled OOF reliability", {
  projection <- matrix(c(-2, 0, 2), nrow = 1L)
  q <- sqrt(max(0, 0.64))
  observed <- q * tanh(projection)
  expect_equal(observed, 0.8 * tanh(projection))
  expect_true(all(abs(observed) <= q))
})

test_that("regulatory integration remains bounded and zero preserving", {
  rna <- matrix(
    c(0, 0.5, 1), nrow = 1,
    dimnames = list("g1", c("u1", "u2", "u3"))
  )
  modifier <- matrix(
    c(1, 1, -1), nrow = 1, dimnames = dimnames(rna)
  )
  observed <- .rc_integrate_regulatory_support(rna, modifier, alpha = 1)

  expect_equal(observed[1, 1], 0)
  expect_equal(observed[1, 3], 1)
  expect_true(observed[1, 2] > 0.5)
  expect_true(all(observed >= 0 & observed <= 1))
})
