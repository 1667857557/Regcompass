test_that("missing reaction expression receives maximum penalty", {
  expression <- matrix(
    c(0, NA, 3, NA),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("R_gene", "R_EX"), c("m1", "m2"))
  )
  roles <- data.frame(
    reaction_id = c("R_gene", "R_EX"),
    role = c("internal", "exchange"),
    role_source = "test",
    stringsAsFactors = FALSE
  )
  result <- rc_compute_multiome_penalty(expression, reaction_roles = roles)
  expect_equal(unname(result$penalty["R_gene", ]), c(1, 1))
  expect_equal(unname(result$penalty["R_EX", ]), c(1, 1))
  expect_true(result$components$penalty_available["R_gene", "m2"])
  expect_true(result$components$missing_expression_imputed_zero[
    "R_gene", "m2"
  ])
  expect_true(result$components$maximum_expression_penalty_flag[
    "R_gene", "m2"
  ])
  expect_true(result$components$observed_zero_expression_flag[
    "R_gene", "m1"
  ])
})

test_that("regulatory alpha is fixed at one with RNA-only fallback for missing modifiers", {
  rna <- matrix(
    c(0, 0.2, 0.8, 1),
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("m1", "m2"))
  )
  modifier <- matrix(
    c(NA, -1, 1, NA),
    nrow = 2,
    dimnames = dimnames(rna)
  )
  observed <- RegCompassR:::.rc_integrate_regulatory_support(
    rna, modifier, alpha = 1
  )
  fallback <- attr(observed, "rna_only_fallback_mask")
  expect_equal(observed[fallback], rna[fallback])
  expect_equal(observed["g1", "m1"], 0)
  expect_equal(observed["g2", "m2"], 1)
  expect_true(all(
    observed[is.finite(observed)] >= 0 &
      observed[is.finite(observed)] <= 1
  ))
  expect_error(
    RegCompassR:::.rc_integrate_regulatory_support(
      rna, modifier, alpha = 0
    ),
    "requires `regulatory_alpha = 1`"
  )
  expect_error(
    RegCompassR:::.rc_integrate_regulatory_support(
      rna, modifier, alpha = 0.5
    ),
    "requires `regulatory_alpha = 1`"
  )
})

test_that("relative score is display-only and composition dependent", {
  penalty <- matrix(
    c(1, 2, 4),
    nrow = 1,
    dimnames = list("target", c("m1", "m2", "m3"))
  )
  score <- rc_compass_score_from_penalty(
    penalty,
    matrix(TRUE, nrow = 1, ncol = 3, dimnames = dimnames(penalty))
  )
  expect_true(isTRUE(attr(score, "composition_dependent")))
  expect_false(isTRUE(attr(score, "effect_size_eligible")))
  expect_true(isTRUE(attr(score, "display_only")))
})

test_that("normalized-unit EB separates depth normalization from shrinkage weight", {
  counts <- matrix(
    c(0, 0, 10, 100, 0, 10),
    nrow = 2,
    dimnames = list(c("low", "high"), c("m1", "m2", "m3"))
  )
  cell_type <- c(m1 = "T", m2 = "T", m3 = "T")
  latent <- RegCompassR:::.rc_latent_metacell_expression(
    counts,
    library_size = c(1e4, 1e5, 1e6),
    mu_min = 0.1,
    cell_type = cell_type
  )
  expect_identical(dimnames(latent$latent_cpm), dimnames(counts))
  expect_true(all(
    latent$zero_class %in% c(
      "observed_positive",
      "observed_zero_continuous_eb"
    )
  ))
  expect_true(all(
    latent$posterior_positive_probability >= 0 &
      latent$posterior_positive_probability <= 1
  ))
  expect_true(all(
    latent$posterior_zero_probability >= 0 &
      latent$posterior_zero_probability <= 1
  ))
  expect_identical(
    latent$model,
    "normalized_unit_gamma_empirical_bayes_by_cell_type_v3"
  )
  expect_identical(latent$prior_estimation_scope, "gene_by_cell_type")
  expect_identical(
    latent$posterior_update_scope,
    "one_normalized_metacell_unit_independent_of_library_size"
  )
  expect_false(latent$hard_structural_zero_used)
  expect_equal(
    latent$observation_weight[, "m1"],
    latent$observation_weight[, "m3"],
    tolerance = 1e-12
  )
  expect_equal(
    latent$prior_weight[, "m1"],
    latent$prior_weight[, "m3"],
    tolerance = 1e-12
  )
})

test_that("equivalent CPM observations receive equivalent EB updates across depth", {
  counts <- matrix(
    c(1, 10, 4, 40),
    nrow = 1,
    dimnames = list("g", c("m1", "m2", "m3", "m4"))
  )
  depth <- c(1e4, 1e5, 2e4, 2e5)
  latent <- RegCompassR:::.rc_latent_metacell_expression(
    counts,
    library_size = depth,
    cell_type = stats::setNames(rep("T", 4), colnames(counts))
  )
  expect_equal(latent$latent_cpm[, "m1"], latent$latent_cpm[, "m2"],
               tolerance = 1e-12)
  expect_equal(latent$latent_cpm[, "m3"], latent$latent_cpm[, "m4"],
               tolerance = 1e-12)
  expect_equal(
    latent$posterior_positive_probability[, "m1"],
    latent$posterior_positive_probability[, "m2"],
    tolerance = 1e-12
  )
})

test_that("zero evidence remains continuous and is never threshold-clamped", {
  counts <- matrix(
    c(0, 0, 1, 3), nrow = 1,
    dimnames = list("g", c("m1", "m2", "m3", "m4"))
  )
  latent <- RegCompassR:::.rc_latent_metacell_expression(
    counts,
    library_size = c(1e4, 1e6, 1e4, 1e4),
    structural_zero_probability = 0.99,
    cell_type = stats::setNames(rep("T", 4), colnames(counts))
  )
  expect_true(all(latent$latent_cpm > 0))
  expect_identical(
    unname(latent$zero_class[, c("m1", "m2")]),
    rep("observed_zero_continuous_eb", 2)
  )
  expect_false(any(latent$zero_class == "credible_structural_zero"))
})

test_that("latent expression rejects normalized non-count input", {
  normalized <- matrix(
    c(0, 0.25, 1.5, 2), nrow = 2,
    dimnames = list(c("g1", "g2"), c("m1", "m2"))
  )
  expect_error(
    RegCompassR:::.rc_latent_metacell_expression(
      normalized,
      library_size = c(1000, 1000),
      cell_type = c(m1 = "T", m2 = "T")
    ),
    "raw integer-like counts"
  )
})

test_that("canonical microCOMPASS defaults keep sample aggregation optional", {
  expect_identical(eval(formals(rc_run_microcompass)$unit)[[1L]], "metacell")
  expect_null(eval(formals(rc_run_microcompass)$sample_col))
  expect_error(
    rc_layer2_unit_matrices(
      list(
        reaction_expression = matrix(
          1, 1, 1, dimnames = list("R1", "m1")
        ),
        unit_meta = data.frame(
          pool_id = "m1", condition = "A", cell_type = "T"
        )
      ),
      unit = "sample_celltype",
      sample_col = NULL,
      celltype_col = "cell_type",
      condition_col = "condition"
    ),
    "requires an explicit"
  )
})

test_that("comparison support policies remain explicit", {
  expect_identical(
    eval(formals(rc_regcompass_step_layer1)$comparison_support),
    c("auto", "pairwise_common", "global_common")
  )
})
