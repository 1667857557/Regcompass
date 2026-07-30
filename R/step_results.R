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
  condition_full_comparison <- .rc_condition_penalty_route(
    layer2, layer2$penalty_condition_full,
    condition_col = params$condition_col,
    celltype_col = params$celltype_col
  )
  rna_only_comparison <- .rc_condition_penalty_route(
    layer2, layer2$penalty_rna_only,
    condition_col = params$condition_col,
    celltype_col = params$celltype_col
  )
  unique_increment_summary <- .rc_condition_increment_summary(
    layer2,
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
  result <- list(
    schema_version = "regcompass_regulatory_metabolic_result_v1",
    version = "2.1.0",
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
    condition_full_exploratory_summary = condition_full_comparison$summary,
    condition_full_exploratory_contrast = condition_full_comparison$contrast,
    rna_only_control_summary = rna_only_comparison$summary,
    rna_only_control_contrast = rna_only_comparison$contrast,
    unique_grn_increment_summary = unique_increment_summary,
    inference_policy = comparison$inference_policy %||%
      "metacell statistical units within one dataset",
    gem_fingerprint = .rc_stage_gem_fingerprint(gem),
    params = list(
      n_conditions = length(conditions),
      requested_condition_col = params$requested_condition_col,
      effective_condition_col = params$condition_col,
      fallback_reason = params$fallback_reason,
      workflow_order = c(
        "single_cell_grn", "native_supercell_metacells", "meta_modules",
        "layer1", "medium_specific_union_gem_layer2"
      ),
      pando_grouping = params$celltype_col,
      pando_design = if (identical(mode, "condition_grn")) {
        paste(
          "shared candidate dictionary, equal-condition transforms, nested",
          "outer-heldout projection, and common support"
        )
      } else {
        "original Pando infer_grn per broad cell type; no condition coefficients"
      },
      pando_regulatory_projection = layer1$projection_provenance,
      metacell_grouping = c(params$condition_col, params$celltype_col),
      metacell_celltype_assignment = "SuperCell cell.annotation",
      metacell_condition_assignment = "SuperCell cell.split.condition",
      metacell_temporary_combined_stratum = FALSE,
      metacell_gamma = params$metacell_args$gamma,
      sample_variable = "not_used",
      meta_module_core_definition =
        "active_pando_targets_complete_gpr_by_effective_group",
      meta_module_expansion =
        "core_subsystem_plus_kegg_reactome_master_rhea_only",
      feasibility_completion = if (
        identical(layer2$model_mode, "meta_module_gem")
      ) {
        "single_global_fastcore_on_each_medium_specific_union_gem"
      } else {
        "not_applicable_full_gem"
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
