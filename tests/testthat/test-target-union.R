target_union_test_gem <- function() {
  S <- diag(5)
  dimnames(S) <- list(paste0("M", 1:5), paste0("R", 1:5))
  reaction_meta <- data.frame(
    reaction_id = paste0("R", 1:5),
    subsystem = c("A", "A", "B", "C", "D"),
    metabolic_module = c("A", "A", "B", "C", "D"),
    kegg_reaction_id = c("K1", NA, "K1", NA, NA),
    reactome_reaction_id = c(NA, NA, NA, NA, NA),
    rhea_master_id = c(NA, NA, "RM1", "RM1", NA),
    stringsAsFactors = FALSE
  )
  gem <- rc_make_gem(
    S,
    lb = rep(0, 5),
    ub = rep(1000, 5),
    reaction_meta = reaction_meta
  )
  gem$gpr_table <- data.frame(
    reaction_id = c("R1", "R1", "R2", "R3", "R4"),
    and_group_id = c("1", "1", "1", "1", "1"),
    gene = c("G1", "G2", "G1", "G3", "G4"),
    stringsAsFactors = FALSE
  )
  gem
}

test_that("selected cores are the only score targets in the expanded union", {
  definition <- .rc_build_target_union_definition(
    gem = target_union_test_gem(),
    core_reaction_ids = "R1"
  )

  expect_identical(definition$global_core_reactions$reaction_id, "R1")
  expect_setequal(
    definition$global_reaction_membership$reaction_id,
    c("R1", "R2", "R3", "R4")
  )
  target <- definition$global_reaction_membership[
    definition$global_reaction_membership$score_target,
    , drop = FALSE
  ]
  support <- definition$global_reaction_membership[
    definition$global_reaction_membership$model_only,
    , drop = FALSE
  ]
  expect_identical(target$reaction_id, "R1")
  expect_setequal(support$reaction_id, c("R2", "R3", "R4"))
  expect_equal(definition$summary$n_score_targets, 1)
  expect_equal(definition$summary$n_model_only_reactions, 3)
})

test_that("gene cores respect complete GPR alternatives by default", {
  gem <- target_union_test_gem()

  one_gene <- .rc_target_union_core_rows(gem, core_genes = "G1")
  expect_identical(one_gene$reaction_id, "R2")
  expect_identical(one_gene$selection_source, "gene_complete_gpr")

  complete <- .rc_target_union_core_rows(
    gem,
    core_genes = c("G1", "G2")
  )
  expect_setequal(complete$reaction_id, c("R1", "R2"))

  direct <- .rc_target_union_core_rows(
    gem,
    core_genes = "G1",
    gene_match = "any_direct"
  )
  expect_setequal(direct$reaction_id, c("R1", "R2"))
})

test_that("invalid manual core selections fail before LP construction", {
  gem <- target_union_test_gem()
  expect_error(
    .rc_target_union_core_rows(gem),
    "Supply at least one"
  )
  expect_error(
    .rc_target_union_core_rows(gem, core_reaction_ids = "missing"),
    "absent from the GEM"
  )
  expect_error(
    .rc_target_union_core_rows(gem, core_genes = "missing"),
    "do not map to GEM GPR rules"
  )
})
