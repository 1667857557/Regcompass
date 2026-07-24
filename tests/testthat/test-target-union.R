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
    S, lb = rep(0, length(reactions)), ub = rep(1000, length(reactions)),
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

target_union_previous_core <- function() {
  data.frame(
    sample_id = "global", module_id = "GLOBAL_UNION",
    reaction_id = c("R1", "R2"), is_core = TRUE,
    stringsAsFactors = FALSE
  )
}

target_union_previous_membership <- function() {
  data.frame(
    sample_id = "global", module_id = "GLOBAL_UNION",
    reaction_id = paste0("R", 1:7),
    is_core = c(TRUE, TRUE, rep(FALSE, 5)),
    inclusion_stage = c(
      "global_union_core", "global_union_core",
      "global_union_biological_member",
      "global_union_biological_member",
      "global_union_biological_member",
      "global_union_biological_member",
      "global_union_local_fastcore_support"
    ),
    stringsAsFactors = FALSE
  )
}

target_union_layer2_stub <- function(file) {
  answer <- list(
    model_mode = "meta_module_gem",
    model_cache_summary = data.frame(
      medium_scenario = "physiologic",
      file = file,
      stringsAsFactors = FALSE
    )
  )
  class(answer) <- c("regcompass_layer2_step", "list")
  answer
}

test_that("only direct database-linked non-core reactions are targets", {
  definition <- .rc_build_target_union_definition(
    gem = target_union_test_gem(),
    global_core_reactions = target_union_previous_core(),
    global_reaction_membership = target_union_previous_membership(),
    core_reaction_ids = "R1"
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
    definition$expanded_scoring_targets$previous_union_is_core
  ))
  core_catalog <- definition$expanded_reaction_catalog[
    definition$expanded_reaction_catalog$previous_union_is_core,
    , drop = FALSE
  ]
  expect_identical(unique(core_catalog$reaction_id), "R2")
  expect_false(any(core_catalog$score_target))
  expect_true(all(
    core_catalog$lp_exclusion_reason ==
      "already_scored_in_original_layer2"
  ))
  expect_equal(definition$summary$n_selected_previous_core, 1)
  expect_equal(definition$summary$n_direct_crossref_relations, 4)
  expect_equal(definition$summary$n_direct_crossref_reactions, 4)
  expect_equal(definition$summary$n_previous_core_reactions_not_rescored, 1)
  expect_equal(definition$summary$n_expanded_score_targets, 3)
  expect_identical(
    definition$summary$expansion_policy,
    "direct_from_selected_core_via_kegg_reactome_master_rhea_only"
  )
  expect_identical(
    definition$summary$scoring_policy,
    "direct_database_crossref_noncore_reactions_only"
  )
  expect_identical(
    definition$summary$model_policy,
    "reuse_exact_previous_global_union_gem"
  )
})

test_that("target-union API has no subsystem or recursive expansion controls", {
  retired <- c("subsystem_table", "expansion_mode", "max_iterations")
  expect_false(any(retired %in% names(formals(rc_regcompass_step_target_union))))
})

test_that("gene selection resolves only previous core anchors", {
  gem <- target_union_test_gem()
  available <- target_union_previous_core()$reaction_id
  one_gene <- .rc_target_union_core_rows(
    gem, available_core_reactions = available, core_genes = "G1"
  )
  expect_identical(one_gene$reaction_id, "R2")
  complete <- .rc_target_union_core_rows(
    gem, available_core_reactions = available,
    core_genes = c("G1", "G2")
  )
  expect_setequal(complete$reaction_id, c("R1", "R2"))
  direct <- .rc_target_union_core_rows(
    gem, available_core_reactions = available,
    core_genes = "G1", gene_match = "any_direct"
  )
  expect_setequal(direct$reaction_id, c("R1", "R2"))
})

test_that("target cache contains only direct non-core database links", {
  gem <- target_union_test_gem()
  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file), add = TRUE)
  saveRDS(gem, file)
  cache <- .rc_build_target_union_model_cache(
    layer2 = target_union_layer2_stub(file),
    target_reactions = c("R3", "R4", "R5"),
    target_direction = "forward"
  )
  expect_length(cache, 3)
  expect_setequal(
    vapply(cache, `[[`, character(1), "reaction_id"),
    c("R3", "R4", "R5")
  )
  expect_true(all(vapply(
    cache, function(x) identical(x$file, file), logical(1)
  )))
  summary <- attr(cache, "summary")
  expect_true(all(summary$reused_without_rebuilding))
  expect_identical(summary$file, file)
  expect_identical(summary$source_model_md5, unname(tools::md5sum(file)))
})

test_that("second LP pass evaluates a non-core target on the cached union", {
  skip_if_not(requireNamespace("highs", quietly = TRUE))
  S <- matrix(
    c(-1, 1, 1, -1), nrow = 2,
    dimnames = list(c("M1", "M2"), c("R1", "R2"))
  )
  gem <- rc_make_gem(
    S,
    lb = c(R1 = 0, R2 = 0),
    ub = c(R1 = 1000, R2 = 1000),
    reaction_meta = data.frame(
      reaction_id = c("R1", "R2"),
      role = "internal", role_source = "test",
      stringsAsFactors = FALSE
    )
  )
  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file), add = TRUE)
  saveRDS(gem, file)
  cache <- .rc_build_target_union_model_cache(
    layer2 = target_union_layer2_stub(file),
    target_reactions = "R2",
    target_direction = "forward"
  )
  layer1 <- list(
    reaction_expression = matrix(
      c(4, 1, 1, 4), nrow = 2, byrow = TRUE,
      dimnames = list(c("R1", "R2"), c("U1", "U2"))
    ),
    unit_meta = data.frame(
      pool_id = c("U1", "U2"),
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
    time_limit = 60,
    flux_threshold = 1e-8,
    parallel = FALSE,
    BPPARAM = FALSE
  )
  expect_equal(dim(result$penalty), c(1, 2))
  expect_true(all(result$evaluated))
  expect_true(all(result$feasible))
  expect_true(all(is.finite(result$penalty)))
  expect_identical(result$target_direction$reaction_id, "R2")
  expect_identical(result$model_cache_summary$file, file)
  expect_true(result$params$structural_model_reused_exactly)
  expect_identical(result$model_file_manifest$file, file)
})

test_that("invalid selections fail before scoring", {
  gem <- target_union_test_gem()
  available <- target_union_previous_core()$reaction_id
  expect_error(
    .rc_target_union_core_rows(gem, available_core_reactions = available),
    "Supply at least one"
  )
  expect_error(
    .rc_target_union_core_rows(
      gem, available_core_reactions = available,
      core_reaction_ids = "missing"
    ),
    "absent from the GEM"
  )
  expect_error(
    .rc_target_union_core_rows(
      gem, available_core_reactions = available,
      core_reaction_ids = "R3"
    ),
    "not core reactions in the previous LP analysis"
  )
  expect_error(
    .rc_target_union_core_rows(
      gem, available_core_reactions = available,
      core_genes = "missing"
    ),
    "do not map to GEM GPR rules"
  )
})
