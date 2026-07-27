#' Assemble compact final RegCompass results
#'
#' The final object contains primary analysis tables, compact regulatory and
#' complete-GPR summaries, scored-reaction annotations, and the compact Layer 2
#' score object required by downstream directional comparisons. Detailed Stage
#' 1--5 intermediates remain in their stage checkpoints and are not duplicated
#' in `regcompass_result.rds`.
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
      !identical(.rc_workflow_signature(grn), .rc_workflow_signature(metacells))) {
    stop("Upstream stages use different workflow parameters.", call. = FALSE)
  }
  .rc_require_stage_gem(grn, gem, "grn")
  .rc_require_stage_gem(meta_modules, gem, "meta_modules")
  .rc_validate_layer1_stage(
    layer1, workflow_params = params, gem = gem, argument = "layer1"
  )
  .rc_validate_layer2_stage(
    layer2, layer1 = layer1, workflow_params = params, gem = gem,
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
  grn_mode <- as.character(
    grn$grn_result$grn_mode %||% grn$params$grn_mode %||%
      "legacy_condition_pando"
  )
  multitask <- identical(grn_mode, "multitask_shared_backbone")
  condition_modules <- meta_modules$condition_modules
  bootstrap_policy <- grn$grn_result$bootstrap_policy %||% list(
    resampling_unit = "cell",
    sample_col = NULL,
    fallback_reason = "legacy Stage 1 object without bootstrap provenance"
  )

  # Full Layer 1 and Layer 2 objects are used transiently to build formal
  # reaction annotations and evidence. Only the compact score subset is retained
  # in the final result; complete stage objects remain in their checkpoint RDS.
  result <- list(
    schema_version = if (multitask) {
      "regcompass_compact_multitask_result_v3"
    } else {
      "regcompass_compact_legacy_result_v2"
    },
    version = "1.8.10",
    species = species,
    model_mode = layer2$model_mode,
    analysis_mode = comparison$analysis_mode,
    grn_mode = grn_mode,
    layer1 = layer1,
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
        "layer1", "medium_specific_union_gem_layer2", "compact_results"
      ),
      grn_mode = grn_mode,
      grn_background = if (multitask) {
        "one_celltype_shared_validated_pando_tf_peak_target_universe"
      } else {
        "independent_condition_x_celltype_pando_candidates"
      },
      grn_condition_model = if (multitask) {
        "global_backbone_plus_symmetric_sum_zero_condition_deviation"
      } else {
        "independent_condition_models"
      },
      grn_stability_policy = if (multitask) {
        paste(
          paste0(
            "condition_stratified_", bootstrap_policy$resampling_unit,
            "_nonparametric_bootstrap;"
          ),
          "stable_projection_weight_equals_full_data_effect_times_selection",
          "frequency_times_conditional_sign_stability"
        )
      } else {
        "legacy_adjusted_p_value_filter"
      },
      sample_col = bootstrap_policy$sample_col %||%
        (grn$params$sample_col %||% NULL),
      bootstrap_resampling_unit =
        bootstrap_policy$resampling_unit %||% NA_character_,
      bootstrap_fallback_reason =
        bootstrap_policy$fallback_reason %||% NA_character_,
      metacell_grouping = params$condition_col,
      metacell_label = params$celltype_col,
      metacell_contract =
        "split_by_condition_then_SuperCell2_label_exact_celltype",
      metacell_gamma = params$metacell_args$gamma,
      condition_balance = if (multitask) {
        "equal_total_GRN_loss_weight_per_condition"
      } else {
        "not_applicable"
      },
      meta_module_core_definition = if (multitask) {
        "condition_celltype_bootstrap_active_targets_complete_gpr"
      } else {
        "condition_celltype_significant_targets_complete_gpr"
      },
      meta_module_expansion =
        "core_subsystem_plus_kegg_reactome_master_rhea_only",
      feasibility_completion = if (identical(layer2$model_mode, "meta_module_gem")) {
        "single_global_fastcore_on_each_medium_specific_union_gem"
      } else {
        "not_applicable_full_gem"
      },
      structural_comparability = paste(
        "all conditions and metacells in one medium reuse identical reaction",
        "IDs, stoichiometric matrix, lower bounds, and upper bounds"
      ),
      penalty_formula = "1/(1+log2(1+E_multiome))",
      result_storage_policy = paste(
        "primary_tables_plus_compact_layer2_scores;",
        "detailed_stage_intermediates_not_embedded"
      ),
      execution_mode = "stepwise"
    )
  )

  result <- .rc_ra_attach_to_result(
    result = result,
    gem = gem,
    condition_col = params$condition_col,
    celltype_col = params$celltype_col
  )

  result$reaction_ranking <- .rc_compact_reaction_ranking(
    result$reaction_ranking
  )
  result$condition_contrast <- .rc_compact_condition_contrast(
    result$condition_contrast
  )
  result$active_regulatory_edges <- .rc_compact_active_edges(
    grn$grn_result, params$condition_col, params$celltype_col
  )
  result$condition_target_genes <- .rc_compact_condition_targets(
    grn$grn_result, params$condition_col, params$celltype_col
  )
  result$core_reactions <- .rc_compact_core_reactions(
    condition_modules, params$condition_col, params$celltype_col
  )
  result$meta_module_summary <- .rc_compact_meta_module_summary(
    condition_modules, params$condition_col, params$celltype_col
  )
  result$grn_metacell_group_coverage <- meta_modules$group_coverage
  result$reaction_catalog <- .rc_compact_reaction_catalog(result$reaction_catalog)
  result$reaction_evidence <- .rc_compact_reaction_evidence(result$reaction_evidence)
  result$microcompass <- .rc_compact_microcompass(layer2)

  result$condition_summary <- NULL
  result$layer1 <- NULL
  result$stage_provenance <- list(
    detailed_intermediates_embedded = FALSE,
    grn_stage_class = class(grn)[[1L]],
    metacell_stage_class = class(metacells)[[1L]],
    meta_module_stage_class = class(meta_modules)[[1L]],
    layer1_stage_class = class(layer1)[[1L]],
    layer2_stage_class = class(layer2)[[1L]],
    compact_layer2_omitted_fields = attr(
      result$microcompass, "omitted_stage5_fields"
    ),
    bootstrap_policy = bootstrap_policy,
    detailed_sources = .rc_result_intermediate_policy()
  )
  result$table_manifest <- .rc_result_table_manifest(list(
    reaction_ranking = result$reaction_ranking,
    condition_contrast = result$condition_contrast,
    active_regulatory_edges = result$active_regulatory_edges,
    condition_target_genes = result$condition_target_genes,
    core_reactions = result$core_reactions,
    meta_module_summary = result$meta_module_summary,
    grn_metacell_group_coverage = result$grn_metacell_group_coverage,
    reaction_catalog = result$reaction_catalog,
    reaction_evidence = result$reaction_evidence
  ))

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  tables <- list(
    reaction_ranking = result$reaction_ranking,
    condition_contrast = result$condition_contrast,
    active_regulatory_edges = result$active_regulatory_edges,
    condition_target_genes = result$condition_target_genes,
    core_reactions = result$core_reactions,
    meta_module_summary = result$meta_module_summary,
    grn_metacell_group_coverage = result$grn_metacell_group_coverage,
    reaction_catalog = result$reaction_catalog,
    reaction_evidence_by_condition_celltype = result$reaction_evidence,
    result_table_manifest = result$table_manifest,
    result_intermediate_policy = result$stage_provenance$detailed_sources
  )
  for (name in names(tables)) {
    value <- tables[[name]]
    if (is.data.frame(value)) {
      .rc_write_tsv_gz(value, file.path(outdir, paste0(name, ".tsv.gz")))
    }
  }

  result <- .rc_step_monitor_finish(result, monitor)
  saveRDS(comparison, file.path(outdir, "step_comparison.rds"))
  saveRDS(result, file.path(outdir, "regcompass_result.rds"))
  result
}
