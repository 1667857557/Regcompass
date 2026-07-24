test_that("target-union scoring admits only annotation-expanded non-core reactions", {
  definition <- .rc_build_target_union_definition(
    gem = target_union_test_gem(),
    global_core_reactions = target_union_previous_core(),
    global_reaction_membership = target_union_previous_membership(),
    core_reaction_ids = "R1"
  )
  allowed_stages <- c(
    "same_core_subsystem",
    "shared_kegg_or_reactome_reaction",
    "shared_master_rhea_reaction"
  )
  expect_true(all(
    definition$expanded_scoring_targets$inclusion_stage %in% allowed_stages
  ))
  expect_false(any(
    definition$expanded_scoring_targets$inclusion_stage %in%
      c("core_grn_gene", "local_fastcore_support", "global_union_local_fastcore_support")
  ))
  expect_false(any(
    definition$expanded_scoring_targets$previous_union_is_core
  ))
  expect_setequal(
    definition$expanded_scoring_targets$reaction_id,
    c("R3", "R4")
  )
})
