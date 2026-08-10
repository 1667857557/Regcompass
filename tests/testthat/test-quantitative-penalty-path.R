test_that("quantitative Pando correction preserves unbounded expression scale", {
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

  expect_equal(as.numeric(out), as.numeric(expected))
  expect_identical(dimnames(out), dimnames(expected))
  expect_equal(out["g1", "u1"], 0)
  expect_gt(out["g2", "u2"], 1)
})

test_that("SuperCell quantitative RNA averages per-cell linear CPM equally", {
  counts <- Matrix::Matrix(
    matrix(
      c(
        10, 100,
        0, 0,
        990, 19900
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("g1", "g2", "other"), c("c1", "c2"))
    ),
    sparse = TRUE
  )
  normalized <- RegCompassR:::.rc_single_cell_linear_cpm(
    counts, genes = c("g1", "g2")
  )
  membership <- data.frame(
    cell_id = c("c1", "c2"),
    metacell_id = c("u1", "u1"),
    stringsAsFactors = FALSE
  )
  averaged <- RegCompassR:::.rc_equal_mean_supercell_expression(
    normalized$expression, membership, units = "u1"
  )

  expect_equal(normalized$expression["g1", "c1"], 10000)
  expect_equal(normalized$expression["g1", "c2"], 5000)
  expect_equal(averaged$expression["g1", "u1"], 7500)
  expect_equal(averaged$expression["g2", "u1"], 0)
  pooled_cpm <- (10 + 100) / (1000 + 20000) * 1e6
  expect_false(isTRUE(all.equal(
    averaged$expression["g1", "u1"], pooled_cpm
  )))
  expect_identical(averaged$library_size_weighted, FALSE)
})

test_that("Pando RNA sources exactly partition the SuperCell membership", {
  validate <- RegCompassR:::.rc_validate_pando_rna_cell_partition

  expect_invisible(validate(
    list(condition = c("c1", "c2"), standard = "c3"),
    c("c3", "c2", "c1")
  ))
  expect_error(
    validate(list(condition = c("c1", "c2")), c("c1")),
    "missing=0, extra=1"
  )
  expect_error(
    validate(list(condition = "c1"), c("c1", "c2")),
    "missing=1, extra=0"
  )
  expect_error(
    validate(list(condition = "c1", standard = "c1"), "c1"),
    "more than one routed Pando RNA source"
  )
})

test_that("Layer 1 v6 quantitative path does not use latent CPM", {
  body_text <- paste(
    deparse(body(RegCompassR:::.rc_cell_first_projection_layer1_v6)),
    collapse = "\n"
  )
  expect_match(body_text, ".rc_quantitative_supercell_rna", fixed = TRUE)
  expect_match(
    body_text,
    "equal_mean_single_cell_linear_cpm",
    fixed = TRUE
  )
  expect_false(grepl(
    "gene_expression_quantitative_rna <- latent$latent_cpm",
    body_text,
    fixed = TRUE
  ))
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

test_that("Layer 2 routes quantitative matrices by explicit marker", {
  structural_multiome <- matrix(
    c(0.2, 0.8), ncol = 1,
    dimnames = list(c("R1", "R2"), "u1")
  )
  structural_rna <- structural_multiome
  attr(
    structural_multiome,
    "regcompass_quantitative_penalty_route"
  ) <- "multiome"
  attr(
    structural_rna,
    "regcompass_quantitative_penalty_route"
  ) <- "rna_only"
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
  expect_identical(
    primary$penalty_route,
    "quantitative_supercell_mean_cpm_multiome"
  )

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
  expect_identical(
    rna_only$penalty_route,
    "quantitative_supercell_mean_cpm_rna_only"
  )

  missing_marker <- layer1
  attr(
    missing_marker$reaction_expression,
    "regcompass_quantitative_penalty_route"
  ) <- NULL
  expect_error(
    RegCompassR:::rc_layer2_unit_matrices(
      missing_marker,
      unit = "metacell",
      sample_col = NULL,
      celltype_col = "cell_type",
      condition_col = "condition"
    ),
    "explicit quantitative penalty route marker"
  )
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
    "compass_quantitative_expression_penalty_v6"
  )
})
