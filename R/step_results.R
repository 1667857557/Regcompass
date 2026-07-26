#' Assemble final RegCompass results
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
      !identical(
        .rc_workflow_signature(grn), .rc_workflow_signature(metacells)
      )) {
    stop("Upstream stages use different workflow parameters.", call. = FALSE)
  }
  .rc_require_stage_gem(grn, gem, "grn")
  .rc_require_stage_gem(meta_modules, gem, "meta_modules")
  .rc_validate_layer1_stage(
    layer1, workflow_params = params, gem = gem, argument = "layer1"
  )
  .rc_validate_layer2_stage(
    layer2,
    layer1 = layer1,
    workflow_params = params,
    gem = gem,
    argument = "layer2"
  )
  species <- .rc_infer_gem_species(gem, species)
  comparison <- .rc_condition_penalty_comparison(
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
  grn_mode <- grn$grn_result$grn_mode %||%
    grn$params$grn_mode %||% "legacy_condition_pando"
  multitask <- identical(grn_mode, "multitask_shared_backbone")
  result <- list(
    schema_version = if (multitask) {
      "regcompass_multitask_shared_grn_v1"
    } else {
      "regcompass_legacy_condition_pando_v1"
    },
    version = "1.8.8",
    species = species,
    model_mode = layer2$model_mode,
    analysis_mode = comparison$analysis_mode,
    grn = grn$grn_result,
    metacells = metacells$pooled,
    layer1 = layer1,
    condition_grn_meta_modules = condition_modules,
    merged_grn_meta_modules = meta_modules$merged_modules,
    grn_meta_modules = meta_modules$merged_modules,
    grn_metacell_group_coverage = meta_modules$group_coverage,
    microcompass = layer2,
    reaction_ranking = comparison$ranking,
    condition_summary = comparison$summary,
    condition_contrast = comparison$contrast,
    inference_policy = comparison$inference_policy,
    gem_fingerprint = .rc_stage_gem_fingerprint(gem),
    params = list(
      n_conditions = length(conditions),
      workflow_order = c(
        "single_cell_grn", "condition_metacells", "meta_modules",
        "layer1", "medium_specific_union_gem_layer2"
      ),
      grn_mode = grn_mode,
      grn_grouping = if (multitask) {
        paste0(
          "one shared candidate universe per ", params$celltype_col,
          "; condition is a joint model task"
        )
      } else {
        c(params$condition_col, params$celltype_col)
      },
      pando_peak_cor =
        grn$grn_result$normalization_policy$pando_peak_cor %||% NA_real_,
      pando_regions = grn$grn_result$normalization_policy$pando_regions,
      metacell_grouping = params$condition_col,
      metacell_celltype_assignment =
        "supercell_label_guided_then_dominant_membership_audit",
      metacell_gamma = params$metacell_args$gamma,
      sample_weighting = if (multitask) {
        "equal total GRN fitting loss per condition"
      } else {
        "none"
      },
      meta_module_core_definition =
        "condition_celltype_active_regulatory_targets_complete_gpr",
      meta_module_expansion =
        "core_subsystem_plus_kegg_reactome_master_rhea_only",
      meta_module_merge = "reaction_id_deduplication_only_not_a_gem",
      feasibility_completion = if (
        identical(layer2$model_mode, "meta_module_gem")
      ) {
        "single_global_fastcore_on_each_medium_specific_union_gem"
      } else {
        "not_applicable_full_gem"
      },
      feasibility_completion_stages = if (
        identical(layer2$model_mode, "meta_module_gem")
      ) {
        "layer2_medium_specific_only"
      } else {
        "none"
      },
      union_gem_definition = paste(
        "medium-constrained merged biological meta-modules plus",
        "global FASTCORE support"
      ),
      structural_comparability = paste(
        "all conditions reuse the exact same medium-specific S, lower bounds,",
        "upper bounds and target catalogue"
      ),
      second_pass_model_policy =
        "reuse_exact_final_medium_specific_union_gem_cache",
      pando_normalization_policy = grn$grn_result$normalization_policy,
      penalty_formula = "1/(1+log2(1+E_multiome))",
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
    result$reaction_catalog,
    file.path(outdir, "reaction_catalog.tsv.gz")
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
