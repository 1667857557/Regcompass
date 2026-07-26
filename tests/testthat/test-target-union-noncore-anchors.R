target_union_noncore_test_gem <- function() {
  reactions <- c("R1", "R2", "R3", "R4")
  S <- diag(length(reactions))
  dimnames(S) <- list(paste0("M", seq_along(reactions)), reactions)
  gem <- rc_make_gem(
    S,
    lb = rep(0, length(reactions)),
    ub = rep(1000, length(reactions)),
    reaction_meta = data.frame(
      reaction_id = reactions,
      subsystem = c("A", "B", "C", "D"),
      metabolic_module = c("A", "B", "C", "D"),
      kegg_reaction_id = c("K1", "K2", "K2", "K1"),
      reactome_reaction_id = NA_character_,
      rhea_master_id = c(NA, "RM2", "RM2", NA),
      role = "internal",
      role_source = "test",
      stringsAsFactors = FALSE
    )
  )
  gem$gpr_table <- data.frame(
    reaction_id = reactions,
    and_group_id = "1",
    gene = c("G1", "G2", "G3", "G4"),
    stringsAsFactors = FALSE
  )
  gem
}

target_union_noncore_core_table <- function() {
  data.frame(
    catalogue_id = "MERGED_META_MODULES",
    reaction_id = "R1",
    is_core = TRUE,
    stringsAsFactors = FALSE
  )
}

target_union_noncore_membership <- function() {
  data.frame(
    catalogue_id = "MERGED_META_MODULES",
    reaction_id = c("R1", "R2", "R3", "R4"),
    is_core = c(TRUE, FALSE, FALSE, FALSE),
    inclusion_stage = c(
      "merged_meta_module_core",
      rep("merged_meta_module_biological_member", 3)
    ),
    stringsAsFactors = FALSE
  )
}

test_that("reaction IDs outside the original core can be remap anchors", {
  definition <- .rc_build_target_union_definition(
    gem = target_union_noncore_test_gem(),
    merged_core_reactions = target_union_noncore_core_table(),
    merged_reaction_membership = target_union_noncore_membership(),
    core_reaction_ids = "R2",
    cached_reaction_ids = c("R1", "R2", "R3", "R4")
  )

  expect_identical(definition$selected_anchor_reactions$reaction_id, "R2")
  expect_false(definition$selected_anchor_reactions$is_core)
  expect_identical(
    definition$selected_anchor_reactions$anchor_role,
    "gem_noncore"
  )
  expect_equal(nrow(definition$selected_core_reactions), 0L)
  expect_identical(definition$selected_noncore_reactions$reaction_id, "R2")
  expect_identical(definition$expanded_scoring_targets$reaction_id, "R3")
  expect_identical(
    unique(definition$expanded_reaction_catalog$anchor_reaction_id),
    "R2"
  )
  expect_false(any(
    definition$expanded_reaction_catalog$anchor_is_original_core
  ))
  expect_equal(definition$summary$n_selected_anchors, 1L)
  expect_equal(definition$summary$n_selected_core, 0L)
  expect_equal(definition$summary$n_selected_noncore_anchors, 1L)
})

test_that("genes can resolve non-core GEM reactions as remap anchors", {
  selected <- .rc_target_union_core_rows(
    gem = target_union_noncore_test_gem(),
    available_core_reactions = "R1",
    core_genes = "G2",
    gene_match = "complete_gpr"
  )

  expect_identical(selected$reaction_id, "R2")
  expect_false(selected$is_core)
  expect_identical(selected$anchor_role, "gem_noncore")
  expect_identical(selected$selection_source, "gene_complete_gpr")
})

test_that("original core anchors retain compatibility outputs", {
  definition <- .rc_build_target_union_definition(
    gem = target_union_noncore_test_gem(),
    merged_core_reactions = target_union_noncore_core_table(),
    merged_reaction_membership = target_union_noncore_membership(),
    core_reaction_ids = "R1",
    cached_reaction_ids = c("R1", "R2", "R3", "R4")
  )

  expect_identical(definition$selected_anchor_reactions$reaction_id, "R1")
  expect_identical(definition$selected_core_reactions$reaction_id, "R1")
  expect_equal(nrow(definition$selected_noncore_reactions), 0L)
  expect_identical(definition$expanded_scoring_targets$reaction_id, "R4")
})
