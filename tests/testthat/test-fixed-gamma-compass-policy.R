test_that("metacell construction uses grouped WNN and gamma 30", {
  builder <- paste(
    deparse(body(RegCompassR:::.rc_build_grouped_wnn_membership)),
    collapse = "\n"
  )
  wrapper <- paste(
    deparse(body(RegCompassR:::.rc_make_condition_celltype_metacells)),
    collapse = "\n"
  )
  defaults <- RegCompassR:::.rc_condition_metacell_defaults()

  expect_identical(defaults$gamma, 30L)
  expect_match(builder, "SCimplify_by_graph_group", fixed = TRUE)
  expect_match(builder, "cell.graph.group", fixed = TRUE)
  expect_match(builder, "cell.split.condition", fixed = TRUE)
  expect_match(
    builder,
    "assay = c\\(rna_assay,[[:space:]]*atac_assay\\)"
  )
  expect_false(grepl("cell.annotation", builder, fixed = TRUE))
  expect_false(grepl("depth_balance", wrapper, fixed = TRUE))
  expect_false(grepl("stratum_col", wrapper, fixed = TRUE))
  expect_false(exists(
    ".rc_native_supercell_membership",
    envir = asNamespace("RegCompassR"),
    inherits = FALSE
  ))
})

test_that("RNA empirical-Bayes priors are estimated by cell type", {
  counts <- matrix(
    c(1, 2, 100, 120, 2, 3, 80, 90),
    nrow = 2,
    byrow = TRUE,
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
    rna, modifier, alpha = 1
  )
  expect_equal(integrated["g1", "u1"], rna["g1", "u1"])
  expect_true(attr(integrated, "rna_only_fallback_mask")["g1", "u1"])
  expect_error(
    RegCompassR:::.rc_integrate_regulatory_support(
      rna, modifier, alpha = 0.25
    ),
    "requires `regulatory_alpha = 1`"
  )
  expect_error(
    rc_regcompass_step_layer1(
      grn = list(), metacells = list(), meta_modules = list(),
      gem = list(), outdir = tempfile(), regulatory_alpha = 0.5
    ),
    "requires `regulatory_alpha = 1`"
  )
})

test_that("GPR aggregation follows COMPASS AND and OR semantics", {
  expect_equal(rc_and_capacity(c(0.2, NA_real_), method = "min"), 0.2)
  expect_equal(rc_and_capacity(c(0.2, NA_real_), method = "mean"), 0.2)
  expect_equal(rc_and_capacity(c(0.2, NA_real_), method = "median"), 0.2)
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
  expect_identical(
    out$missing_expression_policy,
    "compass_missing_expression_max_penalty"
  )
})

test_that("condition contracts are absolute, unversioned and condition-full", {
  expect_identical(
    RegCompassR:::.RC_PANDO_CONDITION_GRN_FIT_SCHEMA,
    "pando_condition_grn_fit"
  )
  validator <- paste(
    deparse(body(RegCompassR:::.rc_require_pando_condition_grn_fit)),
    collapse = "\n"
  )
  extraction <- paste(
    deparse(body(RegCompassR:::.rc_extract_condition_grn_contract)),
    collapse = "\n"
  )
  expect_match(validator, "coefficient_estimable_mask", fixed = TRUE)
  expect_match(validator, "projectable_structural_zero_mask", fixed = TRUE)
  expect_match(validator, "projection_support_mask", fixed = TRUE)
  expect_match(validator, "projection_condition_full_oof", fixed = TRUE)
  expect_match(extraction, "projectable_structural_zero", fixed = TRUE)
  expect_match(extraction, "projection_supported", fixed = TRUE)
  expect_false(grepl("comparison_mask", extraction, fixed = TRUE))
})

test_that("condition Layer 1 records condition-full structural-zero policy", {
  body_text <- paste(
    deparse(body(RegCompassR:::.rc_cell_first_projection_layer1)),
    collapse = "\n"
  )
  expect_match(body_text, "condition_full_oof", fixed = TRUE)
  expect_match(body_text, "gene_projection_condition_unique_oof", fixed = TRUE)
  expect_match(body_text, "nonestimable_projection_policy", fixed = TRUE)
  expect_match(body_text, "condition_grn", fixed = TRUE)
  expect_match(body_text, "standard_pando", fixed = TRUE)
})
