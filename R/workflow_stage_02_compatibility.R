# Compatibility refinements loaded after workflow_stage_01_architecture.R.

# Regulatory evidence follows the same conservative GPR semantics as RNA
# evidence: required subunits use min and alternative isoenzymes use max.
.rc_pando_reaction_confidence <- function(meta_modules, pando_object, gem,
                                           atac_assay = "ATAC",
                                           rna_assay = "RNA") {
  parsed <- rc_parse_gpr_table(gem$gpr_table)
  target_genes <- meta_modules$target_metabolic_genes %||%
    rc_metabolic_gpr_genes(gem$gpr_table)
  gene <- .rc_pando_gene_confidence(
    meta_modules$tf_peak_gene_significant,
    pando_object,
    atac_assay = atac_assay,
    rna_assay = rna_assay,
    target_genes = target_genes
  )
  supported <- rowSums(is.finite(gene$gene_confidence)) > 0L
  gene_for_reaction <- gene$gene_confidence[supported, , drop = FALSE]
  reaction <- rc_reaction_confidence(
    parsed,
    gene_confidence = gene_for_reaction,
    and_method = "min",
    or_method = "max",
    unit_ids = colnames(pando_object)
  )
  list(
    gene_confidence = gene$gene_confidence,
    gene_confidence_diagnostics = gene$diagnostics,
    reaction_confidence = reaction,
    reaction_confidence_matrix = .rc_pando_reaction_confidence_matrix(
      reaction, names(parsed), colnames(pando_object)
    ),
    confidence_source = "pando_signed_tf_peak_gene_regulatory_support"
  )
}

# A pure custom-medium request must not first require a generic exchange
# baseline. Mixed custom and sensitivity requests continue through the canonical
# scenario builder.
rc_make_medium_scenarios <- local({
  canonical <- rc_make_medium_scenarios
  function(
      gem,
      scenario = "permissive_all_exchange",
      custom_medium = NULL,
      uptake_scale = c(
        permissive_all_exchange = 1,
        normal_human_plasma = 1, rpmi1640 = 1, minimal = 0.1,
        low_glucose = 0.1,
        low_glutamine = 0.1, high_lactate = 1
      ),
      condition_col = NULL,
      exchange_roles = c("exchange"),
      condition = condition_col) {
    scenario <- as.character(scenario)
    if (length(scenario) == 1L && identical(scenario, "custom")) {
      if (is.null(custom_medium)) {
        stop("`custom_medium` is required when `scenario = 'custom'`.",
             call. = FALSE)
      }
      required <- c("medium_scenario_id", "exchange_reaction_id",
                    "lb", "ub", "available")
      missing <- setdiff(required, colnames(custom_medium))
      if (length(missing)) {
        stop("`custom_medium` missing columns: ",
             paste(missing, collapse = ", "), call. = FALSE)
      }
      output <- custom_medium
      output$exchange_reaction_id <- trimws(
        as.character(output$exchange_reaction_id)
      )
      output$lb <- suppressWarnings(as.numeric(output$lb))
      output$ub <- suppressWarnings(as.numeric(output$ub))
      output$available <- as.logical(output$available)
      if (anyNA(output$exchange_reaction_id) ||
          any(!nzchar(output$exchange_reaction_id)) ||
          any(!is.finite(output$lb)) || any(!is.finite(output$ub)) ||
          any(output$lb > output$ub) || anyNA(output$available)) {
        stop("Custom medium rows require valid reaction IDs, logical availability and finite ordered bounds.",
             call. = FALSE)
      }
      optional <- c(
        "metabolite_id", "condition", "evidence_source",
        "assumption_level", "target_exchange_flag",
        "concentration_used_for_rate_bound", "rate_bound_source"
      )
      for (name in setdiff(optional, colnames(output))) output[[name]] <- NA
      output$evidence_source[is.na(output$evidence_source)] <-
        "user_supplied_custom_medium"
      output$assumption_level[is.na(output$assumption_level)] <-
        "user_supplied"
      output$concentration_used_for_rate_bound[
        is.na(output$concentration_used_for_rate_bound)
      ] <- FALSE
      output$rate_bound_source[is.na(output$rate_bound_source)] <-
        "user_supplied"
      rownames(output) <- NULL
      return(output)
    }
    canonical(
      gem = gem,
      scenario = scenario,
      custom_medium = custom_medium,
      uptake_scale = uptake_scale,
      condition_col = condition_col,
      exchange_roles = exchange_roles,
      condition = condition
    )
  }
})

# Differential analysis uses raw LP penalty whenever available. This keeps the
# inferential metric in the model's native units instead of testing a display
# rank. Existing score-only objects remain supported.
.rc_original_describe_microcompass_by_group <- rc_describe_microcompass_by_group
rc_describe_microcompass_by_group <- function(
    result,
    sample_col = "sample_id",
    condition_col = "condition",
    celltype_col = "cell_type") {
  metric <- if (!is.null(result$penalty)) "penalty" else "score"
  input <- result
  if (identical(metric, "penalty")) input$score <- result$penalty
  output <- .rc_original_describe_microcompass_by_group(
    input,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col
  )
  output$analysis_metric <- metric
  output$metric_direction <- if (identical(metric, "penalty")) {
    "lower_is_more_supported"
  } else {
    "higher_is_more_supported"
  }
  output
}

.rc_original_test_microcompass_differential <- rc_test_microcompass_differential
rc_test_microcompass_differential <- function(
    result, formula = score ~ condition,
    method = c("lm", "limma_continuous", "wilcox"),
    sample_col = "sample_id", celltype_col = "cell_type",
    condition_col = "condition", covariates = NULL,
    min_samples_per_group = 3, preferred_min_samples_per_group = 5,
    p_adjust_method = "BH", strict_replicate_design = TRUE,
    test_type = c("omnibus", "pairwise")) {
  metric <- if (!is.null(result$penalty)) "penalty" else "score"
  input <- result
  if (identical(metric, "penalty")) input$score <- result$penalty
  output <- .rc_original_test_microcompass_differential(
    input,
    formula = formula,
    method = method,
    sample_col = sample_col,
    celltype_col = celltype_col,
    condition_col = condition_col,
    covariates = covariates,
    min_samples_per_group = min_samples_per_group,
    preferred_min_samples_per_group = preferred_min_samples_per_group,
    p_adjust_method = p_adjust_method,
    strict_replicate_design = strict_replicate_design,
    test_type = test_type
  )
  if (nrow(output)) {
    output$analysis_metric <- metric
    output$metric_direction <- if (identical(metric, "penalty")) {
      "positive_effect_means_higher_penalty_and_weaker_support"
    } else {
      "positive_effect_means_higher_relative_support"
    }
  }
  output
}
