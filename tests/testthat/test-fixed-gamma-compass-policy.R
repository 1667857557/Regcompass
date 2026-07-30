test_that("metacell construction uses one fixed gamma without depth rejection", {
  f <- formals(rc_make_supercell2_metacells)
  expect_identical(eval(f$gamma), 30L)
  expect_identical(eval(f$depth_balance), FALSE)
  body_text <- paste(deparse(body(rc_make_supercell2_metacells)), collapse = "\n")
  expect_match(body_text, "diagnostic_only_no_top1_rejection", fixed = TRUE)
  expect_match(body_text, "global_fixed_across_all_strata", fixed = TRUE)
  expect_false(grepl(".rc_assert_depth_balance", body_text, fixed = TRUE))
})

test_that("RNA empirical-Bayes priors are estimated by cell type", {
  counts <- matrix(
    c(1, 2, 100, 120, 2, 3, 80, 90),
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("a1", "a2", "b1", "b2"))
  )
  depth <- rep(1e6, 4)
  cell_type <- c(a1 = "A", a2 = "A", b1 = "B", b2 = "B")
  out <- RegCompassR:::.rc_latent_metacell_expression(
    counts, depth, cell_type = cell_type
  )
  expect_identical(out$prior_estimation_scope, "gene_by_cell_type")
  expect_identical(colnames(out$prior_mean), c("A", "B"))
  expect_true(all(out$prior_mean[, "B"] > out$prior_mean[, "A"]))
  expect_identical(unname(out$prior_cell_type), unname(cell_type))
})

test_that("regulatory alpha is fixed at one and missing Pando is RNA-only", {
  expect_identical(eval(formals(rc_regcompass_step_layer1)$regulatory_alpha), 1)
  rna <- matrix(
    c(0.2, 0.8), nrow = 1,
    dimnames = list("g1", c("u1", "u2"))
  )
  modifier <- matrix(
    c(NA_real_, 0.5), nrow = 1,
    dimnames = dimnames(rna)
  )
  integrated <- RegCompassR:::.rc_integrate_regulatory_support(
    rna, modifier, alpha = 0.25
  )
  expect_equal(integrated["g1", "u1"], rna["g1", "u1"])
  expect_identical(attr(integrated, "regulatory_alpha"), 1)
  expect_true(attr(integrated, "rna_only_fallback_mask")["g1", "u1"])
  expect_error(
    rc_regcompass_step_layer1(
      grn = list(), metacells = list(), meta_modules = list(),
      gem = list(), outdir = tempfile(), regulatory_alpha = 0.5
    ),
    "requires `regulatory_alpha = 1`"
  )
})

test_that("GPR aggregation follows COMPASS AND and OR semantics", {
  expect_true(is.na(rc_and_capacity(c(0.2, NA_real_), method = "min")))
  expect_equal(rc_and_capacity(c(0.2, NA_real_), method = "mean"), 0.1)
  expect_equal(rc_and_capacity(c(0.2, NA_real_), method = "median"), 0.1)
  expect_equal(rc_or_capacity(c(0.2, NA_real_, 0.7), method = "sum"), 0.9)
  parsed <- list(c("g1", "g2"), "g3")
  score <- c(g1 = 0.4, g2 = 0.2, g3 = 0.7)
  expect_equal(
    rc_reaction_capacity_one(
      parsed, score, and_method = "min", or_method = "sum"
    ),
    0.9
  )
})

test_that("missing reaction expression receives maximum expression penalty", {
  expression <- matrix(
    c(NA_real_, 0, 3), ncol = 1,
    dimnames = list(c("R_missing", "R_zero", "R_high"), "u1")
  )
  out <- rc_compute_multiome_penalty(expression)
  expect_equal(out$penalty["R_missing", "u1"], 1)
  expect_equal(out$penalty["R_zero", "u1"], 1)
  expect_lt(out$penalty["R_high", "u1"], 1)
  expect_true(out$components$penalty_available["R_missing", "u1"])
  expect_true(out$components$missing_expression_imputed_zero[
    "R_missing", "u1"
  ])
  expect_true(out$components$maximum_expression_penalty_flag[
    "R_missing", "u1"
  ])
})

test_that("reference-condition fields are removed from RegCompass contracts", {
  x <- list(
    reference_condition = "Control",
    contrast = 1,
    comparison_mask = TRUE,
    beta_condition_std = matrix(1),
    response_transform = data.frame(
      target = "G1", center = 0, scale = 1,
      reference_condition = "Control"
    )
  )
  out <- RegCompassR:::.rc_strip_reference_contract(x)
  expect_false(any(c(
    "reference_condition", "contrast", "comparison_mask"
  ) %in% names(out)))
  expect_false("reference_condition" %in% colnames(out$response_transform))
  expect_identical(out$coefficient_contract, "absolute_condition_effects_only")
})

test_that("Layer 1 contract admits structural-zero edges to main analysis", {
  body_text <- paste(
    deparse(body(RegCompassR:::.rc_cell_first_projection_layer1)),
    collapse = "\n"
  )
  expect_match(
    body_text, "structural_zero_enters_main_analysis", fixed = TRUE
  )
  expect_match(body_text, "enters_main_analysis = TRUE", fixed = TRUE)
  expect_match(body_text, "structural_zero_contribution <- 0", fixed = TRUE)
})
