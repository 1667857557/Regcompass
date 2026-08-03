.rc_load_condition_modules <- function(meta_modules) {
  value <- meta_modules$condition_modules
  if (is.list(value) && is.data.frame(value$reaction_membership)) {
    return(value)
  }
  reference <- meta_modules$condition_modules_ref %||% value
  if (!is.list(reference) ||
      !identical(
        reference$schema_version,
        "regcompass_external_condition_modules_v1"
      )) {
    stop("The condition meta-module reference is invalid.", call. = FALSE)
  }
  file <- as.character(reference$file %||% "")
  if (length(file) != 1L || !nzchar(file) || !file.exists(file)) {
    stop("The external condition meta-module file is unavailable.",
         call. = FALSE)
  }
  checksum <- unname(tools::md5sum(file))
  if (!identical(checksum, as.character(reference$file_checksum))) {
    stop("The external condition meta-module file failed checksum validation.",
         call. = FALSE)
  }
  readRDS(file)
}

#' Assemble the annotated RegCompass result
#'
#' Supports both condition-aware Pando and automatic standard Pando fallback.
#' A single-condition run returns reaction rankings and summaries with an empty
#' condition contrast rather than manufacturing condition coefficients.
#'
#' @export
rc_regcompass_step_results <- function(
    grn, metacells, meta_modules, layer1, layer2, gem, outdir,
    species = c("auto", "human", "mouse"),
    progress = getOption("RegCompassR.progress", TRUE)) {
  monitor <- .rc_step_monitor_start("results", outdir, progress)
  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)

  materialized_meta_modules <- meta_modules
  materialized_meta_modules$condition_modules <-
    .rc_load_condition_modules(meta_modules)
  meta_modules <- materialized_meta_modules

  .rc_require_stage_class(
    grn, "regcompass_grn_step", "grn", "rc_regcompass_step_grn"
  )
  .rc_require_stage_class(
    metacells, "regcompass_metacell_step", "metacells",
    "rc_regcompass_step_metacells"
  )
  .rc_require_stage_class(
    meta_modules, "regcompass_meta_module_step", "meta_modules",
    "rc_regcompass_step_meta_modules"
  )
  params <- metacells$params
  if (!identical(params, meta_modules$workflow_params) ||
      !identical(.rc_workflow_signature(grn),
                 .rc_workflow_signature(metacells))) {
    stop("Upstream stages use different workflow parameters.",
         call. = FALSE)
  }
  .rc_require_stage_gem(grn, gem, "grn")
  .rc_require_stage_gem(meta_modules, gem, "meta_modules")
  .rc_validate_layer1_stage(
    layer1, workflow_params = params, gem = gem, argument = "layer1"
  )
  .rc_validate_layer2_stage(
    layer2, layer1 = layer1, workflow_params = params,
    gem = gem, argument = "layer2"
  )
  species <- .rc_infer_gem_species(gem, species)
  comparison <- .rc_condition_penalty_comparison(
    layer2,
    condition_col = params$condition_col,
    celltype_col = params$celltype_col
  )
  rna_only_comparison <- .rc_condition_penalty_route(
    layer2,
    layer2$penalty_rna_only,
    condition_col = params$condition_col,
    celltype_col = params$celltype_col
  )
  conditions <- unique(as.character(
    metacells$pooled$metacell_meta[[params$condition_col]]
  ))
  condition_fields <- intersect(c(
    "supported_metabolic_genes", "core_gene_reaction",
    "biological_reaction_membership", "reaction_membership",
    "meta_module_summary", "core_definition", "analysis_group_unit",
    "grn_metacell_group_coverage", "feasibility_completion"
  ), names(meta_modules$condition_modules))
  condition_modules <- meta_modules$condition_modules[condition_fields]
  mode <- params$analysis_mode
  metacell_design <- metacells$pooled$input_design
  is_celltype_union <- identical(layer2$model_mode, "meta_module_gem")
  result <- list(
    schema_version = "regcompass_regulatory_metabolic_result_v3",
    version = "2.4.0",
    species = species,
    model_mode = layer2$model_mode,
    analysis_mode = mode,
    comparison_analysis_mode = comparison$analysis_mode,
    condition_coefficients_calculated = identical(mode, "condition_grn"),
    grn = grn$grn_result,
    metacells = metacells$pooled,
    layer1 = layer1,
    condition_grn_meta_modules = condition_modules,
    merged_grn_meta_modules = meta_modules$merged_modules,
    grn_meta_modules = meta_modules$merged_modules,
    grn_metacell_group_coverage = meta_modules$group_coverage,
    microcompass = layer2,
    reaction_comparison_by_metacell = layer2$comparison_table,
    reaction_ranking = comparison$ranking,
    condition_summary = comparison$summary,
    condition_contrast = comparison$contrast,
    rna_only_control_summary = rna_only_comparison$summary,
    rna_only_control_contrast = rna_only_comparison$contrast,
    inference_policy = comparison$inference_policy %||%
      "metacell statistical units within one dataset",
    gem_fingerprint = .rc_stage_gem_fingerprint(gem),
    params = list(
      n_conditions = length(conditions),
      requested_condition_col = params$requested_condition_col,
      effective_condition_col = params$condition_col,
      fallback_reason = params$fallback_reason,
      workflow_order = c(
        "single_cell_grn",
        "celltype_joint_condition_WNN_metacells",
        "celltype_condition_meta_modules",
        "condition_layer1",
        if (is_celltype_union) {
          "celltype_x_medium_union_gem_layer2"
        } else {
          "shared_full_gem_layer2"
        }
      ),
      pando_grouping = params$celltype_col,
      pando_design = if (identical(mode, "condition_grn")) {
        paste(
          "pooled and condition-specific candidate discovery, exact",
          "TF-peak-target union, and one unscaled fixed-dictionary Gaussian",
          "identity GLM per condition"
        )
      } else {
        "original Pando infer_grn per broad cell type; no condition coefficients"
      },
      pando_regulatory_projection = layer1$projection_provenance,
      primary_penalty = "penalty",
      rna_only_control = "penalty_rna_only",
      nonestimable_edge_policy =
        "coefficient_NA_and_zero_realized_penalty_contribution",
      metacell_purity_grouping = c(params$condition_col, params$celltype_col),
      metacell_graph_grouping = params$celltype_col,
      metacell_supercell_api = metacell_design$native_supercell_api,
      metacell_graph_group_argument = metacell_design$graph_group_argument,
      metacell_condition_argument = metacell_design$condition_argument,
      metacell_graph_method = metacell_design$graph_method,
      metacell_graph_scope = metacell_design$graph_scope,
      metacell_condition_scope = metacell_design$condition_scope,
      metacell_membership_split_timing =
        metacell_design$membership_split_timing,
      metacell_modality_weighting = metacell_design$modality_weighting,
      metacell_temporary_combined_stratum = FALSE,
      metacell_gamma = params$metacell_args$gamma,
      sample_variable = "not_used",
      meta_module_core_definition =
        "active_pando_targets_complete_gpr_by_condition_and_cell_type",
      meta_module_expansion =
        "core_subsystem_plus_kegg_reactome_master_rhea_only",
      meta_module_merge_scope = if (is_celltype_union) {
        "conditions_within_cell_type_only"
      } else {
        "not_used_by_full_gem_mode"
      },
      cross_celltype_meta_module_merge = FALSE,
      structural_scope = layer2$params$structural_scope %||%
        if (is_celltype_union) "cell_type_x_medium" else "full_gem_x_medium",
      shared_across_conditions = layer2$params$shared_across_conditions %||%
        is_celltype_union,
      shared_across_cell_types = layer2$params$shared_across_cell_types %||%
        FALSE,
      union_gem_scope = layer2$params$shared_gem_scope %||%
        if (is_celltype_union) {
          "one_union_gem_per_cell_type_per_medium_shared_within_cell_type"
        } else {
          "not_applicable_full_gem"
        },
      feasibility_completion = if (is_celltype_union) {
        "independent_fastcore_for_each_celltype_x_medium_union_gem"
      } else {
        "not_applicable_full_gem"
      },
      vmax_computation_scope = layer2$params$vmax_computation_scope %||%
        if (is_celltype_union) {
          "celltype_model_x_directional_target_once"
        } else {
          "full_gem_x_directional_target_once"
        },
      pando_normalization_policy = grn$grn_result$normalization_policy,
      penalty_formula = "1/(1+log2(1+E_multiome)); missing E:=0",
      execution_mode = "stepwise"
    )
  )
  result <- .rc_ra_attach_to_result(
    result = result,
    gem = gem,
    condition_col = params$condition_col,
    celltype_col = params$celltype_col
  )
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  .rc_write_tsv_gz(
    result$reaction_catalog, file.path(outdir, "reaction_catalog.tsv.gz")
  )
  .rc_write_tsv_gz(
    result$reaction_evidence,
    file.path(outdir, "reaction_evidence_by_condition_celltype.tsv.gz")
  )
  result <- .rc_step_monitor_finish(result, monitor)
  saveRDS(comparison, file.path(outdir, "step_comparison.rds"))
  saveRDS(result, file.path(outdir, "regcompass_result.rds"))
  result
}
