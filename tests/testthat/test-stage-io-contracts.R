test_that("GRN and metacell groups require bidirectional coverage", {
  grn <- list(condition_fit_status = data.frame(
    condition = c("A", "B"),
    cell_type = c("T", "T"),
    status = c("ok", "ok"),
    n_cells = c(100L, 120L),
    n_active_edges = c(10L, 0L),
    stringsAsFactors = FALSE
  ))
  metacells <- data.frame(
    metacell_id = c("A1", "A2", "B1"),
    condition = c("A", "A", "B"),
    cell_type = "T",
    stringsAsFactors = FALSE
  )
  coverage <- .rc_validate_grn_metacell_group_coverage(
    grn, metacells, "condition", "cell_type"
  )
  expect_true(all(coverage$coverage_complete))
  expect_true(coverage$has_active_pando_evidence[coverage$condition == "A"])
  expect_false(coverage$has_active_pando_evidence[coverage$condition == "B"])
  expect_equal(coverage$n_metacells[coverage$condition == "A"], 2L)
  expect_error(
    .rc_validate_grn_metacell_group_coverage(
      grn,
      metacells[metacells$condition == "A", , drop = FALSE],
      "condition", "cell_type"
    ),
    "do not align"
  )
})

test_that("combined-stratum metacell path is absent", {
  expect_false(exists(
    ".rc_assign_metacell_dominant_celltype",
    inherits = TRUE
  ))
  text <- paste(
    deparse(body(.rc_make_condition_celltype_metacells)),
    collapse = "\n"
  )
  native <- paste(
    deparse(body(.rc_native_supercell_membership)),
    collapse = "\n"
  )
  expect_match(
    native, "SCimplify_by_graph_group_from_embedding", fixed = TRUE
  )
  expect_match(native, "cell.graph.group", fixed = TRUE)
  expect_match(native, "cell.split.condition", fixed = TRUE)
  expect_match(native, ".rc_scale_embedding_block_by_group", fixed = TRUE)
  expect_false(grepl("cell.annotation", native, fixed = TRUE))
  expect_false(grepl(".rc_condition_celltype_pool_col", text, fixed = TRUE))
  expect_false(grepl("stratum_col", text, fixed = TRUE))
})

test_that("metacell stage persists required artifacts", {
  text <- paste(deparse(body(rc_regcompass_step_metacells)), collapse = "\n")
  required <- c(
    "metacell_metadata.tsv.gz",
    "metacell_membership.tsv.gz",
    "metacell_celltype_composition.tsv.gz",
    "metacell_celltype_summary.tsv.gz",
    "merged_metacell_object.rds",
    "step_metacells.rds",
    "SCimplify_by_graph_group_from_embedding",
    "one_independent_graph_per_cell_type",
    "all_conditions_joint_within_cell_type_graph"
  )
  expect_true(all(vapply(required, grepl, logical(1), x = text, fixed = TRUE)))
})

test_that("Stage 3 persists supported genes and core reactions", {
  text <- paste(
    deparse(body(.rc_build_condition_meta_modules)), collapse = "\n"
  )
  required <- c(
    "supported_metabolic_genes.tsv.gz",
    "core_gene_reaction.tsv.gz",
    "meta_module_reactions.tsv.gz",
    "condition_meta_modules.rds"
  )
  expect_true(all(vapply(required, grepl, logical(1), x = text, fixed = TRUE)))
  expect_match(text, ".rc_summarize_supported_metabolic_genes", fixed = TRUE)
  expect_false(grepl("rc_project_metabolic_grn", text, fixed = TRUE))
  expect_false("expansion_mode" %in% names(formals(rc_expand_meta_module_reactions)))
  expect_false("max_iterations" %in% names(formals(rc_expand_meta_module_reactions)))
})

test_that("Layer 1 uses condition-full primary schema", {
  body_text <- paste(
    deparse(body(.rc_cell_first_projection_layer1)), collapse = "\n"
  )
  step_text <- paste(deparse(body(rc_regcompass_step_layer1)), collapse = "\n")
  expect_match(body_text, "regcompass_regulatory_layer1_v2", fixed = TRUE)
  expect_match(body_text, "native_SuperCell_metacell", fixed = TRUE)
  expect_match(body_text, "gene_projection_condition_full_oof", fixed = TRUE)
  expect_match(body_text, "gene_projection_common_oof", fixed = TRUE)
  expect_match(body_text, "gene_projection_condition_unique_oof", fixed = TRUE)
  expect_match(body_text, "reaction_expression_condition_full_oof", fixed = TRUE)
  expect_match(body_text, "and_method = gpr_and_method", fixed = TRUE)
  expect_match(body_text, '"standard_pando"', fixed = TRUE)
  expect_match(body_text, '"condition_grn"', fixed = TRUE)
  expect_false(grepl("reaction_expression_depth_matched_rna", body_text, fixed = TRUE))
  expect_false(grepl("reaction_expression_common_depth_interval_rna", body_text, fixed = TRUE))
  expect_false(grepl("reaction_expression_alpha_sensitivity", body_text, fixed = TRUE))
  expect_false(grepl("reaction_zero_support_sensitivity", body_text, fixed = TRUE))
  expect_false(grepl("reaction_link_saturation_sensitivity", body_text, fixed = TRUE))
  expect_match(step_text, "regcompass_layer1_step", fixed = TRUE)
  expect_match(step_text, "gem_fingerprint", fixed = TRUE)
  expect_match(step_text, "workflow_params", fixed = TRUE)
  expect_true("grn" %in% names(formals(rc_regcompass_step_layer1)))
  expect_identical(
    eval(formals(rc_regcompass_step_layer1)$comparison_support)[[1L]],
    "auto"
  )
  expect_identical(
    eval(formals(rc_regcompass_step_layer1)$gpr_and_method),
    c("min", "median", "mean")
  )
})

test_that("Layer 2 and final results use condition-full primary schema", {
  layer2_text <- paste(deparse(body(rc_regcompass_step_layer2)), collapse = "\n")
  result_text <- paste(deparse(body(rc_regcompass_step_results)), collapse = "\n")
  expect_match(layer2_text, ".rc_validate_layer1_stage", fixed = TRUE)
  expect_match(layer2_text, "regcompass_layer2_step", fixed = TRUE)
  expect_match(layer2_text, "penalty_condition_full_oof", fixed = TRUE)
  expect_match(layer2_text, "penalty_common_oof", fixed = TRUE)
  expect_match(layer2_text, "penalty_condition_unique_increment", fixed = TRUE)
  expect_match(layer2_text, "source_core_reactions", fixed = TRUE)
  expect_false(grepl("penalty_depth_matched_rna", layer2_text, fixed = TRUE))
  expect_false(grepl("penalty_common_depth_interval_rna", layer2_text, fixed = TRUE))
  expect_false(grepl("penalty_alpha_sensitivity", layer2_text, fixed = TRUE))
  expect_match(result_text, ".rc_validate_layer2_stage", fixed = TRUE)
  expect_match(
    result_text,
    "regcompass_regulatory_metabolic_result_v2",
    fixed = TRUE
  )
  expect_match(result_text, 'version = "2.2.0"', fixed = TRUE)
  expect_match(result_text, "condition_grn_meta_modules", fixed = TRUE)
  expect_match(result_text, "merged_grn_meta_modules", fixed = TRUE)
  expect_match(result_text, "supported_metabolic_genes", fixed = TRUE)
  expect_match(result_text, "reaction_catalog", fixed = TRUE)
  expect_match(result_text, "reaction_evidence", fixed = TRUE)
  expect_match(result_text, "metacell_graph_scope", fixed = TRUE)
  expect_match(result_text, "metacell_condition_scope", fixed = TRUE)
  expect_match(result_text, "common_support_component_summary", fixed = TRUE)
  expect_match(result_text, "condition_unique_penalty_increment_summary", fixed = TRUE)
})

test_that("stage validators reject reordered or mismatched units", {
  params <- list(
    condition_col = "condition", celltype_col = "cell_type",
    rna_assay = "RNA", atac_assay = "ATAC",
    analysis_mode = "condition_grn"
  )
  layer1 <- list(
    reaction_expression = matrix(
      1, nrow = 1, ncol = 2,
      dimnames = list("R1", c("U1", "U2"))
    ),
    unit_meta = data.frame(
      pool_id = c("U2", "U1"),
      stringsAsFactors = FALSE
    ),
    workflow_params = params,
    gem_fingerprint = "x",
    analysis_mode = "condition_grn"
  )
  class(layer1) <- c("regcompass_layer1_step", "list")
  expect_error(.rc_validate_layer1_stage(layer1), "not identically aligned")
})
