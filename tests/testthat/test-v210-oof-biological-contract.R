test_that("missing reaction expression remains unavailable", {
  expression <- matrix(
    c(0, NA, 3, NA),
    nrow = 2,
    dimnames = list(c("R_gene", "R_EX"), c("m1", "m2"))
  )
  roles <- data.frame(
    reaction_id = c("R_gene", "R_EX"),
    role = c("internal", "exchange"),
    role_source = "test",
    stringsAsFactors = FALSE
  )
  result <- rc_compute_multiome_penalty(expression, reaction_roles = roles)
  expect_equal(result$penalty["R_gene", "m1"], 1)
  expect_true(is.na(result$penalty["R_gene", "m2"]))
  expect_equal(result$penalty["R_EX", ], c(1, 1))
  expect_false(result$components$penalty_available["R_gene", "m2"])
  expect_true(result$components$penalty_available["R_EX", "m2"])
  expect_true(result$components$observed_zero_expression_flag["R_gene", "m1"])
  expect_false(
    result$components$observed_zero_expression_flag["R_gene", "m2"]
  )
})

test_that("alpha zero is exactly the RNA-only route", {
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
    rna, modifier, alpha = 0
  )
  expect_equal(as.numeric(observed), as.numeric(rna), tolerance = 0)
  expect_identical(dimnames(observed), dimnames(rna))
  expect_equal(observed["g1", "m1"], 0)
})

test_that("bounded regulatory integration preserves zeros and availability", {
  rna <- matrix(
    c(0, 0.25, 0.75, 0.8),
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("m1", "m2"))
  )
  modifier <- matrix(
    c(1, -1, NA, 1),
    nrow = 2,
    dimnames = dimnames(rna)
  )
  observed <- RegCompassR:::.rc_integrate_regulatory_support(
    rna, modifier, alpha = 0.5
  )
  expect_equal(observed["g1", "m1"], 0)
  expect_true(is.na(observed["g1", "m2"]))
  expect_true(all(
    observed[is.finite(observed)] >= 0 &
      observed[is.finite(observed)] <= 1
  ))
  expect_lte(
    (observed["g2", "m2"] / (1 - observed["g2", "m2"])) /
      (rna["g2", "m2"] / (1 - rna["g2", "m2"])),
    sqrt(2) + 1e-12
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

test_that("latent expression separates sampling and structural zeros", {
  counts <- matrix(
    c(0, 0, 10, 20, 0, 1),
    nrow = 2,
    dimnames = list(c("low", "high"), c("m1", "m2", "m3"))
  )
  latent <- RegCompassR:::.rc_latent_metacell_expression(
    counts,
    library_size = c(1e4, 1e5, 1e6),
    mu_min = 0.1
  )
  expect_identical(dimnames(latent$latent_cpm), dimnames(counts))
  expect_true(all(
    latent$zero_class %in% c(
      "observed_positive",
      "sampling_limited_zero",
      "credible_structural_zero"
    )
  ))
  expect_true(all(
    latent$posterior_positive_probability >= 0 &
      latent$posterior_positive_probability <= 1
  ))
  expect_match(latent$marginal_count_model, "negative_binomial")
})

test_that("latent expression rejects normalized non-count input", {
  normalized <- matrix(
    c(0, 0.25, 1.5, 2), nrow = 2,
    dimnames = list(c("g1", "g2"), c("m1", "m2"))
  )
  expect_error(
    RegCompassR:::.rc_latent_metacell_expression(
      normalized, library_size = c(1000, 1000)
    ),
    "raw integer-like counts"
  )
})

test_that("canonical microCOMPASS defaults contain no sample variable", {
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

test_that("diagnostic support policies are absent from the public API", {
  expect_identical(
    eval(formals(rc_regcompass_step_layer1)$comparison_support),
    c("auto", "pairwise_common", "global_common")
  )
})

test_that("common-depth sensitivity never invents a pooled fallback", {
  meta <- data.frame(
    condition = c("A", "A", "B", "B"),
    cell_type = "T"
  )
  observed <- RegCompassR:::.rc_common_depth_interval(
    c(10, 12, 100, 120), meta, "condition", "cell_type"
  )
  expect_false(observed$diagnostics$condition_interval_overlap)
  expect_false(observed$diagnostics$fallback_used)
  expect_false(any(observed$keep))
})
