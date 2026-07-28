test_that("an incomplete alternative branch is not marked as a complete core branch", {
  nodes <- data.frame(
    group_id = c("A", "A"),
    module_id = c("M", "M"),
    gene = c("G1", "G2"),
    stringsAsFactors = FALSE
  )
  gpr <- data.frame(
    reaction_id = c("R", "R", "R"),
    and_group_id = c("isozyme_1", "complex_2", "complex_2"),
    gene = c("G1", "G2", "G3"),
    stringsAsFactors = FALSE
  )

  mapped <- RegCompassR:::rc_map_meta_module_core_reactions(nodes, gpr)
  complete <- mapped[mapped$and_group_id == "isozyme_1", , drop = FALSE]
  incomplete <- mapped[mapped$and_group_id == "complex_2", , drop = FALSE]

  expect_true(all(complete$reaction_is_core))
  expect_true(all(complete$group_complete))
  expect_true(all(complete$is_core))
  expect_identical(unique(complete$inclusion_stage), "core_complete_gpr")

  expect_true(all(incomplete$reaction_is_core))
  expect_false(any(incomplete$group_complete))
  expect_false(any(incomplete$is_core))
  expect_true(all(incomplete$is_partial_candidate))
  expect_identical(
    unique(incomplete$inclusion_stage),
    "alternative_incomplete_gpr_branch"
  )

  hard <- RegCompassR:::.rc_hard_core_rows(mapped)
  expect_setequal(hard$gene, "G1")
  expect_setequal(hard$reaction_id, "R")
})
