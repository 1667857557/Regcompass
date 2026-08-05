test_that("native Layer 1 evidence uses condition median then condition max", {
  units <- c(
    "ctrl_1", "ctrl_2", "ctrl_3",
    "treat_1", "treat_2", "treat_3"
  )
  rna <- matrix(
    c(
      1, 1, 100, 2, 2, 2,
      4, 6, 8, 1, 1, 50
    ),
    nrow = 2L, byrow = TRUE,
    dimnames = list(c("R1", "R2"), units)
  )
  multiome <- matrix(
    c(
      3, 3, 99, 4, 4, 4,
      2, 2, 2, 9, 9, 100
    ),
    nrow = 2L, byrow = TRUE,
    dimnames = dimnames(rna)
  )
  regulatory <- matrix(
    c(
      0.1, 0.1, 1, 0.4, 0.4, 0.4,
      0.2, 0.2, 0.2, 0.8, 0.8, 1
    ),
    nrow = 2L, byrow = TRUE,
    dimnames = dimnames(rna)
  )
  layer1 <- list(
    reaction_expression_rna_only = rna,
    reaction_expression = multiome,
    reaction_regulatory_support_fraction = regulatory,
    unit_meta = data.frame(
      unit_id = units,
      condition = rep(c("Control", "Treatment"), each = 3L),
      cell_type = rep("epithelial", length(units)),
      stringsAsFactors = FALSE
    ),
    workflow_params = list(
      condition_col = "condition",
      celltype_col = "cell_type"
    )
  )

  tables <- RegCompassR:::.rc_corda_layer1_evidence_tables(layer1)

  expect_identical(tables$source, "native_layer1_metacell_matrices")
  expect_identical(tables$aggregation, "condition_median_max")
  expect_equal(
    tables$support$rna_reaction_support[
      match(c("R1", "R2"), tables$support$reaction_id)
    ],
    c(2, 6)
  )
  expect_equal(
    tables$multiome$multiome_reaction_support[
      match(c("R1", "R2"), tables$multiome$reaction_id)
    ],
    c(4, 9)
  )
  expect_equal(
    tables$regulatory$regulatory_support[
      match(c("R1", "R2"), tables$regulatory$reaction_id)
    ],
    c(0.4, 0.8)
  )
  expect_true(all(tables$support$n_conditions == 2L))
  expect_true(all(tables$support$n_metacells == 6L))
})

test_that("condition median prevents a single metacell outlier from driving CORDA2", {
  value <- matrix(
    c(1, 1, 100, 2, 2, 2),
    nrow = 1L,
    dimnames = list("R1", paste0("u", 1:6))
  )
  result <- RegCompassR:::.rc_corda_condition_median_max(
    value = value,
    unit_meta = data.frame(unit_id = colnames(value)),
    condition = rep(c("Control", "Treatment"), each = 3L),
    cell_type = rep("T", 6L)
  )

  expect_equal(result$value, 2)
  expect_false(result$value == 100)
})

test_that("legacy CORDA2 reaction-support tables remain supported", {
  layer1 <- list(
    reaction_support = data.frame(
      reaction_id = c("R1", "R2"),
      cell_type = c("T", "T"),
      rna_reaction_support = c(0.2, 0.8),
      multiome_reaction_support = c(0.3, 0.9),
      stringsAsFactors = FALSE
    )
  )

  tables <- RegCompassR:::.rc_corda_layer1_evidence_tables(layer1)

  expect_identical(tables$source, "legacy_reaction_support_tables")
  expect_identical(tables$aggregation, "legacy_celltype_mean")
  expect_identical(tables$support, layer1$reaction_support)
  expect_identical(tables$multiome, layer1$reaction_support)
})

test_that("native CORDA2 evidence rejects misaligned metacell matrices", {
  rna <- matrix(
    1:4,
    nrow = 2L,
    dimnames = list(c("R1", "R2"), c("u1", "u2"))
  )
  multiome <- rna[, c("u2", "u1"), drop = FALSE]
  layer1 <- list(
    reaction_expression_rna_only = rna,
    reaction_expression = multiome,
    unit_meta = data.frame(
      unit_id = c("u1", "u2"),
      condition = c("A", "B"),
      cell_type = c("T", "T"),
      stringsAsFactors = FALSE
    ),
    workflow_params = list(
      condition_col = "condition",
      celltype_col = "cell_type"
    )
  )

  expect_error(
    RegCompassR:::.rc_corda_layer1_evidence_tables(layer1),
    "must align exactly"
  )
})
