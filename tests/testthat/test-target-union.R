target_union_test_gem <- function() {
  reactions <- paste0("R", 1:7)
  S <- diag(length(reactions))
  dimnames(S) <- list(paste0("M", seq_along(reactions)), reactions)
  reaction_meta <- data.frame(
    reaction_id = reactions,
    subsystem = c("A", "B", "C", "D", "E", "A", "F"),
    metabolic_module = c("A", "B", "C", "D", "E", "A", "F"),
    kegg_reaction_id = c("K1", "K1", "K1", NA, NA, NA, NA),
    reactome_reaction_id = c("RE1", NA, NA, "RE1", NA, NA, NA),
    rhea_master_id = c("RM1", NA, "RM2", NA, "RM1", NA, "RM2"),
    role = "internal",
    role_source = "test",
    stringsAsFactors = FALSE
  )
  gem <- rc_make_gem(
    S,
    lb = rep(0, length(reactions)),
    ub = rep(1000, length(reactions)),
    reaction_meta = reaction_meta
  )
  gem$gpr_table <- data.frame(
    reaction_id = c("R1", "R1", "R2", "R3", "R4", "R5", "R6", "R7"),
    and_group_id = "1",
    gene = c("G1", "G2", "G1", "G3", "G4", "G5", "G6", "G7"),
    stringsAsFactors = FALSE
  )
  gem
}

target_union_merged_core <- function() {
  data.frame(
    catalogue_id = "MERGED_META_MODULES",
    reaction_id = c("R1", "R2"),
    is_core = TRUE,
    stringsAsFactors = FALSE
  )
}

target_union_merged_membership <- function() {
  data.frame(
    catalogue_id = "MERGED_META_MODULES",
    reaction_id = paste0("R", 1:7),
    is_core = c(TRUE, TRUE, rep(FALSE, 5)),
    inclusion_stage = c(
      "merged_meta_module_core",
      "merged_meta_module_core",
      rep("merged_meta_module_biological_member", 5)
    ),
    stringsAsFactors = FALSE
  )
}

target_union_layer2_stub <- function(gem, scenario = "physiologic") {
  file <- tempfile(fileext = ".rds")
  gem$is_union_gem <- TRUE
  gem$union_gem_medium_scenario <- scenario
  gem$union_gem_scope <-
    "one_medium_shared_across_conditions_and_metacells"
  gem$target_status <- "ok"
  gem$closure_diagnostics <- data.frame()
  saveRDS(gem, file)
  checksum <- unname(tools::md5sum(file))
  answer <- list(
    model_mode = "meta_module_gem",
    model_cache_summary = data.frame(
      medium_scenario = scenario,
      file = file,
      file_checksum = checksum,
      build_strategy = "medium_specific_union_gem",
      completion_stage = "single_global_fastcore_after_meta_module_merge",
      stringsAsFactors = FALSE
    )
  )
  class(answer) <- c("regcompass_layer2_step", "list")
  list(layer2 = answer, file = file, checksum = checksum)
}

test_that("only direct database-linked non-core reactions are targets", {
  definition <- .rc_build_target_union_definition(
    gem = target_union_test_gem(),
    merged_core_reactions = target_union_merged_core(),
    merged_reaction_membership = target_union_merged_membership(),
    core_reaction_ids = "R1",
    cached_reaction_ids = paste0("R", 1:7)
  )
  expect_identical(definition$selected_core_reactions$reaction_id, "R1")
  expect_setequal(
    definition$expanded_reaction_catalog$reaction_id,
    c("R2", "R3", "R4", "R5")
  )
  expect_setequal(
    definition$expanded_scoring_targets$reaction_id,
    c("R3", "R4", "R5")
  )
  expect_false(any(
    definition$expanded_reaction_catalog$reaction_id %in% c("R6", "R7")
  ))
  expect_setequal(
    definition$expanded_reaction_catalog$expansion_type,
    c(
      "shared_kegg_reaction",
      "shared_reactome_reaction",
      "shared_master_rhea_reaction"
    )
  )
  expect_true(all(definition$expanded_scoring_targets$score_target))
  expect_false(any(
    definition$expanded_scoring_targets$merged_catalogue_is_core
  ))
  core_catalogue <- definition$expanded_reaction_catalog[
    definition$expanded_reaction_catalog$merged_catalogue_is_core,
    , drop = FALSE
  ]
  expect_identical(unique(core_catalogue$reaction_id), "R2")
  expect_false(any(core_catalogue$score_target))
  expect_true(all(
    core_catalogue$lp_exclusion_reason ==
      "already_scored_in_original_layer2"
  ))
  expect_equal(definition$summary$n_selected_core, 1)
  expect_equal(definition$summary$n_direct_crossref_relations, 4)
  expect_equal(definition$summary$n_direct_crossref_reactions, 4)
  expect_equal(definition$summary$n_merged_core_reactions_not_rescored, 1)
  expect_equal(definition$summary$n_expanded_score_targets, 3)
  expect_identical(
    definition$summary$expansion_policy,
    "direct_from_selected_core_via_kegg_reactome_master_rhea_only"
  )
  expect_identical(
    definition$summary$scoring_policy,
    "direct_noncore_reactions_present_in_all_final_union_gems"
  )
  expect_identical(
    definition$summary$model_policy,
    "reuse_exact_final_medium_specific_union_gems"
  )
})

test_that("target-union API has no reconstruction or scoring timeout controls", {
  retired <- c(
    "subsystem_table", "expansion_mode", "max_iterations",
    "fastcore_epsilon", "max_support_reactions", "strict", "time_limit"
  )
  expect_false(any(retired %in% names(formals(rc_regcompass_step_target_union))))
  expect_false("time_limit" %in% names(formals(.rc_score_existing_union_cache)))
  expect_error(
    rc_regcompass_step_target_union(
      layer1 = list(), meta_modules = list(), layer2 = list(), gem = list(),
      outdir = tempfile(), layer2_args = list(time_limit = 60)
    ),
    "output of `rc_regcompass_step_meta_modules"
  )
})

test_that("gene selection resolves only original Layer 2 core anchors", {
  gem <- target_union_test_gem()
  available <- target_union_merged_core()$reaction_id
  one_gene <- .rc_target_union_core_rows(
    gem, available_core_reactions = available, core_genes = "G1"
  )
  expect_identical(one_gene$reaction_id, "R2")
  complete <- .rc_target_union_core_rows(
    gem,
    available_core_reactions = available,
    core_genes = c("G1", "G2")
  )
  expect_setequal(complete$reaction_id, c("R1", "R2"))
  direct <- .rc_target_union_core_rows(
    gem,
    available_core_reactions = available,
    core_genes = "G1",
    gene_match = "any_direct"
  )
  expect_setequal(direct$reaction_id, c("R1", "R2"))
  expect_error(
    .rc_target_union_core_rows(
      gem,
      available_core_reactions = available,
      core_reaction_ids = "R1",
      core_genes = "G5"
    ),
    "do not resolve to original Layer 2 core targets"
  )
})

test_that("final union-GEM membership is authoritative over Stage 3 catalogue", {
  membership <- target_union_merged_membership()
  membership <- membership[membership$reaction_id != "R5", , drop = FALSE]
  definition <- .rc_build_target_union_definition(
    gem = target_union_test_gem(),
    merged_core_reactions = target_union_merged_core(),
    merged_reaction_membership = membership,
    core_reaction_ids = "R1",
    cached_reaction_ids = paste0("R", 1:7)
  )
  expect_true("R5" %in% definition$expanded_scoring_targets$reaction_id)
  r5 <- definition$expanded_reaction_catalog[
    definition$expanded_reaction_catalog$reaction_id == "R5",
    , drop = FALSE
  ]
  expect_true(all(r5$available_in_all_cached_union_gems))
  expect_false(any(r5$present_in_merged_catalogue))
  expect_identical(
    unique(r5$merged_catalogue_inclusion_stage),
    "global_fastcore_support_not_in_merged_catalogue"
  )
})

test_that("target cache reuses exact final union GEM files", {
  stub <- target_union_layer2_stub(target_union_test_gem())
  on.exit(unlink(stub$file), add = TRUE)
  before <- unname(tools::md5sum(stub$file))

  expect_setequal(
    .rc_target_union_cached_reaction_ids(stub$layer2),
    paste0("R", 1:7)
  )
  cache <- .rc_build_target_union_model_cache(
    layer2 = stub$layer2,
    target_reactions = c("R3", "R4", "R5"),
    target_direction = "forward"
  )
  expect_length(cache, 3)
  expect_setequal(
    vapply(cache, `[[`, character(1), "reaction_id"),
    c("R3", "R4", "R5")
  )
  expect_true(all(vapply(
    cache, function(x) identical(x$file, stub$file), logical(1)
  )))
  summary <- attr(cache, "summary")
  expect_true(all(summary$reused_without_rebuilding))
  expect_identical(summary$file, stub$file)
  expect_identical(summary$file_checksum, stub$checksum)
  expect_identical(
    summary$second_pass_build_strategy,
    "reuse_exact_final_medium_specific_union_gem"
  )
  expect_identical(unname(tools::md5sum(stub$file)), before)
})

test_that("tampered or non-union cache files are rejected", {
  stub <- target_union_layer2_stub(target_union_test_gem())
  on.exit(unlink(stub$file), add = TRUE)
  bad_checksum <- stub$layer2
  bad_checksum$model_cache_summary$file_checksum <- "not-the-checksum"
  expect_error(
    .rc_target_union_model_summary(bad_checksum),
    "checksum"
  )

  plain_file <- tempfile(fileext = ".rds")
  on.exit(unlink(plain_file), add = TRUE)
  saveRDS(target_union_test_gem(), plain_file)
  plain <- stub$layer2
  plain$model_cache_summary$file <- plain_file
  plain$model_cache_summary$file_checksum <- unname(tools::md5sum(plain_file))
  expect_error(
    .rc_target_union_model_summary(plain),
    "original final medium-specific union GEM"
  )
})

test_that("second LP pass evaluates targets without rebuilding or timeout", {
  skip_if_not(requireNamespace("highs", quietly = TRUE))
  S <- matrix(
    c(-1, 1, 1, -1),
    nrow = 2,
    dimnames = list(c("M1", "M2"), c("R1", "R2"))
  )
  gem <- rc_make_gem(
    S,
    lb = c(R1 = 0, R2 = 0),
    ub = c(R1 = 1000, R2 = 1000),
    reaction_meta = data.frame(
      reaction_id = c("R1", "R2"),
      role = "internal",
      role_source = "test",
      stringsAsFactors = FALSE
    )
  )
  stub <- target_union_layer2_stub(gem)
  on.exit(unlink(stub$file), add = TRUE)
  before <- unname(tools::md5sum(stub$file))
  cache <- .rc_build_target_union_model_cache(
    layer2 = stub$layer2,
    target_reactions = "R2",
    target_direction = "forward"
  )
  layer1 <- list(
    reaction_expression = matrix(
      c(4, 1, 1, 4),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("R1", "R2"), c("U1", "U2"))
    ),
    unit_meta = data.frame(
      pool_id = c("U1", "U2"),
      unit_id = c("U1", "U2"),
      sample_id = "S1",
      condition = c("A", "B"),
      cell_type = "C",
      stringsAsFactors = FALSE
    )
  )
  result <- .rc_score_existing_union_cache(
    layer1 = layer1,
    gem = gem,
    model_cache = cache,
    condition_col = "condition",
    sample_col = "sample_id",
    celltype_col = "cell_type",
    omega = 0.95,
    solver = "highs",
    flux_threshold = 1e-8,
    parallel = FALSE,
    BPPARAM = FALSE
  )
  expect_equal(dim(result$penalty), c(1, 2))
  expect_true(all(result$evaluated))
  expect_true(all(result$feasible))
  expect_true(all(is.finite(result$penalty)))
  expect_identical(result$target_direction$reaction_id, "R2")
  expect_identical(result$model_cache_summary$file, stub$file)
  expect_identical(result$model_mode, "reused_medium_specific_union_gem")
  expect_true(result$params$structural_model_reused_exactly)
  expect_false(result$params$fastcore_rerun)
  expect_false(result$params$model_rebuild)
  expect_identical(result$params$scoring_time_limit, "none")
  expect_identical(result$model_file_manifest$file, stub$file)
  expect_identical(unname(tools::md5sum(stub$file)), before)
})

test_that("invalid selections fail before scoring", {
  gem <- target_union_test_gem()
  available <- target_union_merged_core()$reaction_id
  expect_error(
    .rc_target_union_core_rows(gem, available_core_reactions = available),
    "Supply at least one"
  )
  expect_error(
    .rc_target_union_core_rows(
      gem,
      available_core_reactions = available,
      core_reaction_ids = "missing"
    ),
    "absent from the GEM"
  )
  noncore <- .rc_target_union_core_rows(
    gem,
    available_core_reactions = available,
    core_reaction_ids = "R3"
  )
  expect_identical(noncore$reaction_id, "R3")
  expect_false(noncore$is_core)
  expect_identical(noncore$anchor_role, "gem_noncore")
  expect_error(
    .rc_target_union_core_rows(
      gem,
      available_core_reactions = available,
      core_genes = "missing"
    ),
    "do not map to GEM GPR rules"
  )
})
