test_that("meta-module expansion is one fixed ordered pass", {
  S <- diag(10)
  dimnames(S) <- list(paste0("M", 1:10), paste0("R", 1:10))
  reaction_meta <- data.frame(
    reaction_id = paste0("R", 1:10),
    subsystem = rep(LETTERS[1:5], each = 2),
    metabolic_module = rep(LETTERS[1:5], each = 2),
    kegg_reaction_id = c("K1", NA, "K1", NA, NA, NA, "K2", NA, "K2", NA),
    reactome_reaction_id = c(NA, "X1", NA, NA, "X1", NA, NA, NA, NA, NA),
    rhea_master_id = c("RM1", NA, NA, NA, "RM2", NA, "RM2", NA, NA, NA),
    stringsAsFactors = FALSE
  )
  gem <- rc_make_gem(
    S,
    lb = rep(0, 10),
    ub = rep(1000, 10),
    reaction_meta = reaction_meta
  )
  core <- data.frame(
    sample_id = "S1",
    module_id = "S1::SUPPORTED_METABOLIC_GENES",
    gene = "G1",
    reaction_id = "R1",
    stringsAsFactors = FALSE
  )

  expanded <- rc_expand_meta_module_reactions(gem, core)

  expect_setequal(
    expanded$reaction_membership$reaction_id,
    c("R1", "R2", "R3", "R5", "R7")
  )
  expect_false("R9" %in% expanded$reaction_membership$reaction_id)
  expect_false(any(
    c("R4", "R6", "R8", "R10") %in%
      expanded$reaction_membership$reaction_id
  ))
  stage <- stats::setNames(
    expanded$reaction_membership$inclusion_stage,
    expanded$reaction_membership$reaction_id
  )
  expect_identical(stage[["R2"]], "same_core_subsystem")
  expect_identical(stage[["R3"]], "shared_kegg_or_reactome_reaction")
  expect_identical(stage[["R7"]], "shared_master_rhea_reaction")
  expect_equal(expanded$summary$n_core_reactions, 1)
  expect_equal(expanded$summary$n_subsystem_added, 1)
  expect_equal(expanded$summary$n_database_added, 2)
  expect_equal(expanded$summary$n_rhea_added, 1)
  expect_identical(
    expanded$summary$expansion_policy,
    "single_ordered_annotation_pass"
  )
  expect_false("expansion_mode" %in% colnames(expanded$reaction_membership))
  expect_false("iterations" %in% colnames(expanded$summary))
  expect_false("expansion_mode" %in% names(formals(rc_expand_meta_module_reactions)))
  expect_false("max_iterations" %in% names(formals(rc_expand_meta_module_reactions)))
  expect_false(exists(".rc_meta_module_one_hop", inherits = TRUE))
})

test_that("partial GPR anchors are removed before module summarization", {
  S <- diag(3)
  dimnames(S) <- list(paste0("M", 1:3), paste0("R", 1:3))
  gem <- rc_make_gem(
    S,
    reaction_meta = data.frame(
      reaction_id = paste0("R", 1:3),
      subsystem = rep("A", 3),
      metabolic_module = rep("A", 3),
      stringsAsFactors = FALSE
    )
  )
  core <- data.frame(
    sample_id = "S1",
    module_id = "S1::SUPPORTED_METABOLIC_GENES",
    gene = c("G1", "G2"),
    reaction_id = c("R1", "R2"),
    and_group_id = "1",
    required_genes = c("G1", "G2;G3"),
    matched_genes = c("G1", "G2"),
    missing_genes = c("", "G3"),
    group_complete = c(TRUE, FALSE),
    is_core = c(TRUE, FALSE),
    is_partial_candidate = c(FALSE, TRUE),
    inclusion_stage = c("core_complete_gpr", "partial_gpr_candidate"),
    stringsAsFactors = FALSE
  )

  expanded <- rc_expand_meta_module_reactions(gem, core)
  expect_setequal(expanded$reaction_membership$reaction_id, c("R1", "R3"))
  expect_equal(expanded$summary$n_core_genes, 1)
  expect_equal(expanded$summary$n_core_reactions, 1)
  expect_equal(expanded$summary$n_reactions, 2)
})

test_that("UNASSIGNED subsystem labels are not pooled", {
  S <- diag(3)
  dimnames(S) <- list(paste0("M", 1:3), paste0("R", 1:3))
  gem <- rc_make_gem(
    S,
    reaction_meta = data.frame(
      reaction_id = paste0("R", 1:3),
      subsystem = c("A", "UNASSIGNED", "UNASSIGNED"),
      metabolic_module = c("A", "UNASSIGNED", "UNASSIGNED"),
      stringsAsFactors = FALSE
    )
  )
  maps <- rc_reaction_crossref_maps(gem)
  expect_false(any(toupper(maps$subsystem$subsystem_id) == "UNASSIGNED"))
})

test_that("active Pando targets define one supported gene set per group", {
  significant <- data.frame(
    group_id = rep("control|tumor", 4),
    condition = rep("control", 4),
    cell_type = rep("tumor", 4),
    tf = c("TF1", "TF2", "TF3", "TF4"),
    target = c("G1", "G1", "G2", "G3"),
    region = paste0("chr1-", 1:4, "-", 11:14),
    estimate = c(1, -0.5, -1, 0.8),
    padj = c(0.001, 0.01, 0.02, 0.03),
    rsq = c(0.4, 0.4, 0.3, 0.2),
    stringsAsFactors = FALSE
  )
  grn_result <- list(
    group_cols = c("condition", "cell_type"),
    tf_peak_gene_condition = significant
  )
  supported <- .rc_summarize_supported_metabolic_genes(
    grn_result,
    metabolic_genes = c("G1", "G2", "G3", "G4")
  )
  expect_setequal(supported$gene, c("G1", "G2", "G3"))
  expect_equal(length(unique(supported$module_id)), 1L)
  expect_equal(
    supported$n_significant_edges[supported$gene == "G1"],
    2L
  )
  expect_equal(
    supported$n_active_edges[supported$gene == "G1"],
    2L
  )
  expect_identical(
    unique(supported$evidence_definition),
    "active_pando_tf_peak_gene_target"
  )
  expect_equal(
    supported$n_positive_edges[supported$gene == "G1"],
    1L
  )
  expect_equal(
    supported$n_negative_edges[supported$gene == "G1"],
    1L
  )

  gpr <- data.frame(
    reaction_id = c("R1", "R1", "R2", "R3", "R3"),
    and_group_id = c(1, 1, 1, 1, 1),
    gene = c("G1", "G2", "G3", "G1", "G4"),
    stringsAsFactors = FALSE
  )
  mapped <- rc_map_meta_module_core_reactions(
    supported[, c("sample_id", "module_id", "gene")],
    gpr
  )
  expect_setequal(
    unique(mapped$reaction_id[mapped$is_core]),
    c("R1", "R2")
  )
  expect_false(unique(mapped$is_core[mapped$reaction_id == "R3"]))
})

test_that("GRN mapping remains explicit from catalogue to final union GEM", {
  skip_if_not(
    requireNamespace("highs", quietly = TRUE) ||
      requireNamespace("Rglpk", quietly = TRUE) ||
      requireNamespace("gurobi", quietly = TRUE)
  )
  solver <- if (requireNamespace("highs", quietly = TRUE)) {
    "highs"
  } else if (requireNamespace("Rglpk", quietly = TRUE)) {
    "glpk"
  } else {
    "gurobi"
  }
  S <- diag(3)
  dimnames(S) <- list(paste0("M", 1:3), paste0("R", 1:3))
  gpr <- data.frame(
    reaction_id = c("R1", "R2", "R3"),
    and_group_id = 1L,
    gene = c("G1", "G2", "G3"),
    stringsAsFactors = FALSE
  )
  gem <- rc_make_gem(
    S,
    lb = rep(0, 3),
    ub = rep(1000, 3),
    reaction_meta = data.frame(
      reaction_id = paste0("R", 1:3),
      metabolic_module = c("A", "A", "B"),
      role = "internal",
      role_source = "curated",
      stringsAsFactors = FALSE
    )
  )
  gem$gpr_table <- gpr
  supported <- data.frame(
    group_id = "S1",
    sample_id = "S1",
    gene = c("G1", "G2"),
    module_id = "S1::SUPPORTED_METABOLIC_GENES",
    stringsAsFactors = FALSE
  )
  core <- rc_map_meta_module_core_reactions(supported, gpr)
  membership <- data.frame(
    sample_id = "S1",
    module_id = "S1::SUPPORTED_METABOLIC_GENES",
    reaction_id = c("R1", "R2"),
    is_core = TRUE,
    stringsAsFactors = FALSE
  )
  condition_modules <- list(
    condition_fit_status = data.frame(group_id = "S1"),
    tf_peak_gene_condition_all = data.frame(),
    tf_peak_gene_condition = data.frame(),
    supported_metabolic_genes = supported,
    core_gene_reaction = core,
    reaction_membership = membership,
    meta_module_summary = data.frame()
  )
  merged <- .rc_merge_meta_module_catalogue(condition_modules)
  expect_false(merged$is_gem)
  expect_setequal(
    merged$merged_reaction_membership$reaction_id,
    c("R1", "R2")
  )

  union_gem <- .rc_complete_medium_union_gem(
    gem = gem,
    reaction_membership = merged$merged_reaction_membership,
    core_reactions = merged$merged_core_reactions,
    solver = solver,
    strict = TRUE
  )
  expect_true(union_gem$is_union_gem)
  expect_setequal(colnames(union_gem$S), c("R1", "R2"))
  expect_true(all(union_gem$reaction_meta$merged_meta_module_member))
  expect_false(any(union_gem$reaction_meta$global_fastcore_support))
})

test_that("Pando validation enforces the configured repository", {
  description <- list(
    RemoteUsername = "1667857557",
    RemoteRepo = "Pando_regcompass",
    RemoteRef = "HEAD",
    RemoteSha = "any-sha"
  )
  valid <- .rc_validate_pando_repository(
    description = description,
    installed_version = "9.9.9"
  )
  expect_identical(valid$version, "9.9.9")
  expect_identical(valid$remote_username, "1667857557")
  expect_identical(valid$remote_repo, "Pando_regcompass")

  bad_user <- description
  bad_user$RemoteUsername <- "other-user"
  expect_error(
    .rc_validate_pando_repository(
      description = bad_user,
      installed_version = "9.9.9"
    ),
    "remote username mismatch"
  )
  bad_repo <- description
  bad_repo$RemoteRepo <- "Pando"
  expect_error(
    .rc_validate_pando_repository(
      description = bad_repo,
      installed_version = "9.9.9"
    ),
    "remote repository mismatch"
  )
})
