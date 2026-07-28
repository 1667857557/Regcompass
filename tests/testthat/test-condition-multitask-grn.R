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

test_that("RegCompass enforces the shared-design independent Pando contract", {
  defaults <- eval(
    formals(.rc_run_condition_single_cell_grns)$pando_infer_args
  )
  expect_identical(defaults$method, "shared_design_independent")
  expect_identical(defaults$candidate_screen, "condition_union")
  expect_identical(defaults$condition_mix, 1)
  expect_identical(defaults$condition_weight, "equal")
  expect_true(defaults$scale)

  body_text <- paste(
    deparse(body(.rc_run_condition_single_cell_grns)), collapse = "\n"
  )
  expect_match(body_text, ".rc_extract_condition_grn_contract", fixed = TRUE)
  expect_match(body_text, "condition_grn_fit_v2.rds", fixed = TRUE)
  expect_match(body_text, "pando_edge_predictor_transforms", fixed = TRUE)
})

test_that("TF-by-ATAC activity uses the Pando interaction predictor", {
  tf <- matrix(
    c(1, 2, 3, 4), nrow = 2,
    dimnames = list(c("TF1", "TF2"), c("u1", "u2"))
  )
  peak <- matrix(
    c(5, 6, 7, 8), nrow = 2,
    dimnames = list(c("P1", "P2"), c("u1", "u2"))
  )

  expect_equal(.rc_tf_peak_interaction(tf, peak), tf * peak)
  expect_error(
    .rc_tf_peak_interaction(tf, peak[, 1, drop = FALSE]),
    "identical dimensions"
  )
})

test_that("model-space projection preserves coefficient-effect magnitude", {
  edge_model <- matrix(
    c(1, -1, 2, 3), nrow = 2, byrow = TRUE,
    dimnames = list(c("edge1", "edge2"), c("u1", "u2"))
  )
  delta <- c(0.5, -0.25)
  observed <- .rc_project_condition_edges(
    edge_model, delta, target_rsq = 0.64
  )
  expected <- sqrt(0.64) * as.numeric(crossprod(delta, edge_model))

  expect_equal(observed, expected)
  expect_equal(
    .rc_project_condition_edges(edge_model, 2 * delta, 0.64),
    2 * observed
  )
  expect_equal(
    .rc_project_condition_edges(edge_model, c(0, 0), 0.64),
    c(0, 0)
  )
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
