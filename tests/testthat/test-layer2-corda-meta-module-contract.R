test_that("CORDA2 resolves Stage 3 merged reaction membership", {
  meta_modules <- list(
    workflow_params = list(celltype_col = "broad_type"),
    merged_modules = list(
      celltype_col = "broad_type",
      merged_reaction_membership = data.frame(
        broad_type = c("epithelial", "epithelial", "immune"),
        reaction_id = c("R1", "R2", "R3"),
        inclusion_stage = c("core", "expanded", "core"),
        stringsAsFactors = FALSE
      )
    )
  )

  membership <- RegCompassR:::.rc_meta_module_reaction_membership(
    meta_modules
  )

  expect_equal(
    membership$cell_type,
    c("epithelial", "epithelial", "immune")
  )
  expect_equal(membership$reaction_id, c("R1", "R2", "R3"))
  expect_true("inclusion_stage" %in% colnames(membership))
})

test_that("CORDA2 evidence mapping executes the complete Stage 3 handoff", {
  units <- c("ctrl_1", "ctrl_2", "treat_1", "treat_2")
  rna <- matrix(
    c(
      1, 3, 7, 9,
      2, 4, 1, 1
    ),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("R1", "R2"), units)
  )
  multiome <- matrix(
    c(
      2, 4, 8, 10,
      3, 5, 2, 2
    ),
    nrow = 2L,
    byrow = TRUE,
    dimnames = dimnames(rna)
  )
  layer1 <- list(
    reaction_expression_rna_only = rna,
    reaction_expression = multiome,
    unit_meta = data.frame(
      unit_id = units,
      condition = c("Control", "Control", "Treatment", "Treatment"),
      broad_type = rep("epithelial", length(units)),
      stringsAsFactors = FALSE
    ),
    workflow_params = list(
      condition_col = "condition",
      celltype_col = "broad_type"
    )
  )
  meta_modules <- list(
    workflow_params = list(celltype_col = "broad_type"),
    merged_modules = list(
      celltype_col = "broad_type",
      merged_reaction_membership = data.frame(
        broad_type = "epithelial",
        reaction_id = "R1",
        stringsAsFactors = FALSE
      )
    )
  )

  evidence <- RegCompassR:::.rc_corda_reaction_evidence(
    layer1 = layer1,
    meta_modules = meta_modules,
    regulatory_weight = 0.20
  )

  r1 <- evidence[evidence$reaction_id == "R1", , drop = FALSE]
  r2 <- evidence[evidence$reaction_id == "R2", , drop = FALSE]
  expect_equal(nrow(r1), 1L)
  expect_equal(nrow(r2), 1L)
  expect_true(r1$merged_meta_module_member)
  expect_false(r2$merged_meta_module_member)
  expect_identical(r1$evidence_aggregation, "condition_median_max")
  expect_equal(r1$rna_support, 8)
  expect_equal(r1$multiome_support, 9)
})

test_that("CORDA2 rejects missing Stage 3 membership explicitly", {
  expect_error(
    RegCompassR:::.rc_meta_module_reaction_membership(list()),
    "does not contain non-empty merged reaction membership"
  )
})
