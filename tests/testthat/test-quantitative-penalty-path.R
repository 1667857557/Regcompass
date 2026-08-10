test_that("quantitative Pando correction preserves latent-CPM scale", {
  rna <- matrix(
    c(0, 1, 10, 100),
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("u1", "u2"))
  )
  modifier <- matrix(
    c(-1, 1, 0.5, NA_real_),
    nrow = 2,
    dimnames = dimnames(rna)
  )

  out <- RegCompassR:::.rc_integrate_regulatory_expression(
    rna, modifier
  )
  expected_modifier <- modifier
  expected_modifier[!is.finite(expected_modifier)] <- 0
  expected <- rna * 2^expected_modifier

  expect_equal(out, expected, ignore_attr = TRUE)
  expect_equal(out["g1", "u1"], 0)
  expect_gt(out["g2", "u2"], 1)
  expect_identical(
    attr(out, "integration_formula"),
    "X_multiome=X_RNA*2^R; X_RNA=latent_CPM; nonfinite R:=0; R clipped to [-1,1]"
  )
})

test_that("bounded support remains separate from quantitative expression", {
  latent_cpm <- matrix(
    c(1, 10, 100, 1000),
    nrow = 1,
    dimnames = list("g1", paste0("u", 1:4))
  )
  bounded <- rc_gene_score(
    log1p(latent_cpm), mode = "absolute", half_saturation = 1
  )

  expect_true(all(bounded >= 0 & bounded < 1))
  expect_gt(max(latent_cpm), 1)

  bounded_penalty <- rc_compute_multiome_penalty(
    bounded,
    reaction_roles = data.frame(
      reaction_id = "g1", role = "internal", stringsAsFactors = FALSE
    )
  )$penalty[1, 4]
  quantitative_penalty <- rc_compute_multiome_penalty(
    latent_cpm,
    reaction_roles = data.frame(
      reaction_id = "g1", role = "internal", stringsAsFactors = FALSE
    )
  )$penalty[1, 4]

  expect_gt(bounded_penalty, 0.5)
  expect_lt(quantitative_penalty, 0.1)
})

test_that("Layer 2 LP prefers quantitative matrices and preserves RNA control", {
  structural_multiome <- matrix(
    c(0.2, 0.8), ncol = 1,
    dimnames = list(c("R1", "R2"), "u1")
  )
  structural_rna <- matrix(
    c(0.1, 0.7), ncol = 1,
    dimnames = dimnames(structural_multiome)
  )
  quantitative_multiome <- matrix(
    c(20, 200), ncol = 1,
    dimnames = dimnames(structural_multiome)
  )
  quantitative_rna <- matrix(
    c(10, 100), ncol = 1,
    dimnames = dimnames(structural_multiome)
  )
  layer1 <- list(
    reaction_expression = structural_multiome,
    reaction_expression_rna_only = structural_rna,
    reaction_expression_quantitative = quantitative_multiome,
    reaction_expression_quantitative_rna_only = quantitative_rna,
    unit_meta = data.frame(
      pool_id = "u1", cell_type = "T", condition = "A",
      stringsAsFactors = FALSE
    )
  )

  primary <- RegCompassR:::rc_layer2_unit_matrices(
    layer1,
    unit = "metacell",
    sample_col = NULL,
    celltype_col = "cell_type",
    condition_col = "condition"
  )
  expect_equal(primary$reaction_expression, quantitative_multiome)
  expect_match(primary$summary, "quantitative_latent_cpm_multiome")

  control <- layer1
  control$reaction_expression <- structural_rna
  rna_only <- RegCompassR:::rc_layer2_unit_matrices(
    control,
    unit = "metacell",
    sample_col = NULL,
    celltype_col = "cell_type",
    condition_col = "condition"
  )
  expect_equal(rna_only$reaction_expression, quantitative_rna)
  expect_match(rna_only$summary, "quantitative_latent_cpm_rna_only")
})

test_that("quantitative COMPASS cost keeps high-expression dynamic range", {
  expression <- matrix(
    c(1, 10, 100, 1000),
    nrow = 1,
    dimnames = list("R1", paste0("u", 1:4))
  )
  out <- rc_compute_multiome_penalty(
    expression,
    reaction_roles = data.frame(
      reaction_id = "R1", role = "internal", stringsAsFactors = FALSE
    )
  )

  expected <- 1 / (1 + log2(1 + expression))
  expect_equal(out$penalty, expected)
  expect_equal(out$penalty[1, 1], 0.5)
  expect_lt(out$penalty[1, 4], 0.1)
  expect_identical(
    out$penalty_version,
    "compass_quantitative_expression_penalty_v5"
  )
})
