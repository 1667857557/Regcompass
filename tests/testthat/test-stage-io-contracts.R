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

canonical_metacell_stage <- function() {
  params <- list(
    condition_col = "condition", celltype_col = "cell_type",
    rna_assay = "RNA", atac_assay = "ATAC"
  )
  design <- list(
    native_supercell_api = "SCimplify_by_graph_group",
    graph_group_argument = "cell.graph.group",
    condition_argument = "cell.split.condition",
    graph_method = "multimodal_WNN",
    graph_scope = "one_independent_WNN_graph_per_cell_type",
    condition_scope = "all_conditions_joint_within_cell_type_graph",
    membership_split_timing = "after_joint_WNN_graph_clustering",
    modality_weighting = "adaptive_WNN_within_cell_type",
    temporary_combined_stratum = FALSE
  )
  contract <- c(list(
    schema_version = "regcompass_celltype_wnn_condition_joint_cache",
    condition_col = params$condition_col,
    celltype_col = params$celltype_col,
    rna_assay = params$rna_assay,
    atac_assay = params$atac_assay
  ), design[setdiff(names(design), c("graph_method", "temporary_combined_stratum"))])
  out <- list(
    params = params,
    pooled = list(input_design = design, cache_contract = contract)
  )
  class(out) <- c("regcompass_metacell_step", "list")
  out
}

test_that("Stage 2 accepts only the canonical grouped-WNN contract", {
  stage <- canonical_metacell_stage()
  expect_invisible(.rc_validate_metacell_artifact_contract(stage))
  old <- stage
  old$pooled$input_design$native_supercell_api <-
    "SCimplify_by_graph_group_from_embedding"
  expect_error(
    .rc_validate_metacell_artifact_contract(old),
    "joint-condition WNN"
  )
})

test_that("Stage 2 uses one grouped multimodal SuperCell API", {
  require_text <- paste(deparse(body(.rc_require_supercell_api)), collapse = "\n")
  build_text <- paste(deparse(body(.rc_build_grouped_wnn_membership)), collapse = "\n")
  expect_match(require_text, "SCimplify_by_graph_group", fixed = TRUE)
  expect_match(build_text, "cell.graph.group", fixed = TRUE)
  expect_match(build_text, "cell.split.condition", fixed = TRUE)
  expect_false(exists(".rc_native_supercell_membership", inherits = TRUE))
  expect_false(exists(".rc_scale_embedding_block_by_group", inherits = TRUE))
})

test_that("metacell stage persists required artifacts", {
  text <- paste(c(
    deparse(body(rc_regcompass_step_metacells)),
    deparse(body(RegCompassR:::.rc_original_step_metacells_cell_set_contract))
  ), collapse = "\n")
  required <- c(
    "metacell_metadata.tsv.gz",
    "metacell_membership.tsv.gz",
    "metacell_celltype_composition.tsv.gz",
    "metacell_celltype_summary.tsv.gz",
    "merged_metacell_object.rds",
    "step_metacells.rds"
  )
  expect_true(all(vapply(required, grepl, logical(1), x = text, fixed = TRUE)))
})

test_that("Layer 1 producer and validator use the same schema", {
  producer <- paste(deparse(body(.rc_cell_first_projection_layer1)), collapse = "\n")
  validator <- paste(deparse(body(.rc_validate_layer1_stage)), collapse = "\n")
  expect_match(producer, "regcompass_regulatory_layer1_v3", fixed = TRUE)
  expect_match(validator, "regcompass_regulatory_layer1_v3", fixed = TRUE)
  expect_false(grepl("regcompass_regulatory_layer1_v2", validator, fixed = TRUE))
})

test_that("Layer 2 and results retain condition-full primary routing", {
  layer2_text <- paste(deparse(body(rc_regcompass_step_layer2)), collapse = "\n")
  result_text <- paste(deparse(body(rc_regcompass_step_results)), collapse = "\n")
  expect_match(layer2_text, ".rc_validate_layer1_stage", fixed = TRUE)
  expect_match(layer2_text, "penalty_condition_full_oof", fixed = TRUE)
  expect_match(result_text, ".rc_validate_layer2_stage", fixed = TRUE)
  expect_match(result_text, 'version = "2.2.4"', fixed = TRUE)
  expect_match(result_text, "metacell_modality_weighting", fixed = TRUE)
})

test_that("stage validators reject reordered units", {
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
    unit_meta = data.frame(pool_id = c("U2", "U1")),
    workflow_params = params,
    gem_fingerprint = "x",
    analysis_mode = "condition_grn"
  )
  class(layer1) <- c("regcompass_layer1_step", "list")
  expect_error(.rc_validate_layer1_stage(layer1), "not identically aligned")
})
