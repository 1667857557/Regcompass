.hierarchy_metacell_artifact <- function() {
  structure(list(
    params = list(
      condition_col = "condition",
      celltype_col = "cell_type",
      rna_assay = "RNA",
      atac_assay = "ATAC"
    ),
    pooled = list(
      partition_policy = "hierarchy_constrained",
      partition_schema_version = "shared_walktrap_condition_cut_v1",
      input_design = list(
        native_supercell_api = "SCimplify_by_graph_group",
        graph_group_argument = "cell.graph.group",
        condition_argument = "cell.split.condition",
        condition_partition = "hierarchy_constrained",
        partition_schema_version = "shared_walktrap_condition_cut_v1",
        graph_method = "multimodal_WNN",
        clustering_method = "one_shared_walktrap_hierarchy_per_cell_type",
        final_partition_method =
          "condition_specific_finest_feasible_cut_of_shared_hierarchy",
        aggregation_method = "SCimplify_for_Seurat_with_membership",
        graph_scope = "one_independent_WNN_graph_per_cell_type",
        condition_scope =
          "all_conditions_joint_for_WNN_and_Walktrap_then_condition_specific_hierarchy_cut",
        membership_split_timing =
          "during_final_shared_hierarchy_cut_selection",
        modality_weighting = "adaptive_WNN_within_cell_type",
        temporary_combined_stratum = FALSE
      ),
      cache_contract = list(
        schema_version = "regcompass_shared_walktrap_condition_cut_cache_v1",
        condition_col = "condition",
        celltype_col = "cell_type",
        rna_assay = "RNA",
        atac_assay = "ATAC",
        native_supercell_api = "SCimplify_by_graph_group",
        graph_group_argument = "cell.graph.group",
        condition_argument = "cell.split.condition",
        condition_partition = "hierarchy_constrained",
        graph_scope = "one_independent_WNN_graph_per_cell_type",
        condition_scope =
          "shared_WNN_and_Walktrap_with_condition_specific_hierarchy_cut",
        membership_split_timing =
          "condition_specific_cut_of_shared_walktrap_hierarchy",
        graph_method = "SuperCell_multimodal_WNN_then_walktrap",
        modality_weighting = "adaptive_WNN_within_cell_type",
        aggregation_method = "SCimplify_for_Seurat_membership_mode"
      )
    )
  ), class = c("regcompass_metacell_step", "list"))
}

test_that("current shared hierarchy Stage 2 artifacts pass validation", {
  artifact <- .hierarchy_metacell_artifact()
  expect_invisible(
    RegCompassR:::.rc_validate_metacell_artifact_contract(artifact)
  )
  expect_invisible(
    RegCompassR:::.rc_require_stage_class(
      artifact, "regcompass_metacell_step", "metacells",
      "rc_regcompass_step_metacells"
    )
  )
})

test_that("legacy post-split Stage 2 contracts are rejected", {
  artifact <- .hierarchy_metacell_artifact()
  artifact$pooled$cache_contract$schema_version <-
    "regcompass_celltype_wnn_condition_joint_cache"
  expect_error(
    RegCompassR:::.rc_validate_metacell_artifact_contract(artifact),
    "shared-Walktrap hierarchy-constrained"
  )
})
