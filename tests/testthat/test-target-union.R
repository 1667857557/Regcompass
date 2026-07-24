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
    role = rep("internal", 5),
    role_source = rep("test", 5),
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

target_union_previous_core <- function() {
  data.frame(
    sample_id = "global",
    module_id = "GLOBAL_UNION",
    reaction_id = c("R1", "R2"),
    is_core = TRUE,
    stringsAsFactors = FALSE
  )
}

target_union_previous_membership <- function() {
  data.frame(
    sample_id = "global",
    module_id = "GLOBAL_UNION",
    reaction_id = paste0("R", 1:5),
    is_core = c(TRUE, TRUE, FALSE, FALSE, FALSE),
    inclusion_stage = c(
      "global_union_core", "global_union_core",
      "global_union_biological_member",
      "global_union_biological_member",
      "global_union_local_fastcore_support"
    ),
    stringsAsFactors = FALSE
  )
}

test_that("selected core and annotation-related reactions are all LP targets", {
  definition <- .rc_build_target_union_definition(
    gem = target_union_test_gem(),
    global_core_reactions = target_union_previous_core(),
    global_reaction_membership = target_union_previous_membership(),
    core_reaction_ids = "R1"
  )

  expect_identical(definition$selected_core_reactions$reaction_id, "R1")
  expect_setequal(
    definition$expanded_scoring_targets$reaction_id,
    c("R1", "R2", "R3", "R4")
  )
  expect_true(all(definition$expanded_scoring_targets$score_target))
  expect_false("model_only" %in% colnames(
    definition$expanded_scoring_targets
  ))
  expect_identical(
    definition$expanded_scoring_targets$reaction_id[
      definition$expanded_scoring_targets$selected_core_anchor
    ],
    "R1"
  )
  expect_equal(definition$summary$n_selected_previous_core, 1)
  expect_equal(definition$summary$n_expanded_score_targets, 4)
  expect_identical(
    definition$summary$model_policy,
    "reuse_previous_global_union_gem_without_rebuilding"
  )
})

test_that("gene selection resolves only previously scored core reactions", {
  gem <- target_union_test_gem()
  available <- target_union_previous_core()$reaction_id

  one_gene <- .rc_target_union_core_rows(
    gem,
    available_core_reactions = available,
    core_genes = "G1"
  )
  expect_identical(one_gene$reaction_id, "R2")
  expect_identical(
    one_gene$selection_source,
    "previous_core_gene_complete_gpr"
  )

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
})

test_that("previous union GEM files are reused for expanded targets", {
  gem <- target_union_test_gem()
  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file), add = TRUE)
  saveRDS(gem, file)
  previous_layer2 <- list(
    model_mode = "meta_module_gem",
    model_cache_summary = data.frame(
      medium_scenario = "physiologic",
      file = file,
      stringsAsFactors = FALSE
    )
  )

  cache <- .rc_build_target_union_model_cache(
    layer2 = previous_layer2,
    target_reactions = c("R1", "R2", "R3", "R4"),
    target_direction = "forward"
  )

  expect_length(cache, 4)
  expect_true(all(vapply(
    cache,
    function(x) identical(x$file, file),
    logical(1)
  )))
  expect_true(all(
    attr(cache, "summary")$reused_without_rebuilding
  ))
  expect_setequal(
    attr(cache, "direction_diagnostics")$reaction_id,
    c("R1", "R2", "R3", "R4")
  )
})

test_that("second LP pass runs against the cached union GEM", {
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
      role = c("internal", "internal"),
      role_source = c("test", "test"),
      stringsAsFactors = FALSE
    )
  )
  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file), add = TRUE)
  saveRDS(gem, file)
  previous_layer2 <- list(
    model_mode = "meta_module_gem",
    model_cache_summary = data.frame(
      medium_scenario = "physiologic",
      file = file,
      stringsAsFactors = FALSE
    )
  )
  cache <- .rc_build_target_union_model_cache(
    layer2 = previous_layer2,
    target_reactions = c("R1", "R2"),
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
      sample_id = c("S1", "S1"),
      condition = c("A", "B"),
      cell_type = c("C", "C"),
      stringsAsFactors = FALSE
    )
  )
  result <- .rc_score_target_union_cache(
    layer1 = layer1,
    gem = gem,
    model_cache = cache,
    medium_scenarios = data.frame(
      medium_scenario_id = "physiologic",
      exchange_reaction_id = NA_character_,
      lb = NA_real_,
      ub = NA_real_,
      available = FALSE,
      .no_constraints = TRUE,
      stringsAsFactors = FALSE
    ),
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

  expect_equal(dim(result$penalty), c(2, 2))
  expect_true(all(result$evaluated))
  expect_true(all(result$feasible))
  expect_true(all(is.finite(result$penalty)))
  expect_setequal(result$target_direction$reaction_id, c("R1", "R2"))
  expect_identical(result$model_mode, "reused_global_union_gem")
})

test_that("invalid selections fail before the second LP pass", {
  gem <- target_union_test_gem()
  available <- target_union_previous_core()$reaction_id
  expect_error(
    .rc_target_union_core_rows(
      gem,
      available_core_reactions = available
    ),
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
  expect_error(
    .rc_target_union_core_rows(
      gem,
      available_core_reactions = available,
      core_reaction_ids = "R3"
    ),
    "not core reactions in the previous LP analysis"
  )
  expect_error(
    .rc_target_union_core_rows(
      gem,
      available_core_reactions = available,
      core_genes = "missing"
    ),
    "do not map to GEM GPR rules"
  )
})
