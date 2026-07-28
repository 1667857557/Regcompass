# Enrich the canonical Stage 1 output contract after inference.

.rc_run_celltype_multitask_grns_core <- .rc_run_celltype_multitask_grns

.rc_run_celltype_multitask_grns <- function(...) {
  call_args <- list(...)
  answer <- do.call(.rc_run_celltype_multitask_grns_core, call_args)
  condition_all <- answer$tf_peak_gene_condition_all
  if (is.data.frame(condition_all) && nrow(condition_all)) {
    condition_col <- answer$group_cols[[1L]]
    celltype_col <- answer$group_cols[[2L]]
    stability_columns <- intersect(c(
      "group_id", condition_col, celltype_col, "edge_universe_id",
      "model_edge_universe_id", "edge_id", "tf", "region",
      "atac_feature_id", "target", "n_observable_conditions",
      "observable_conditions", "global_estimate", "condition_deviation",
      "effective_estimate", "coefficient_parameterization",
      "theta_penalty_factor", "selection_frequency",
      "selection_frequency_mc_se", "selection_frequency_lower_95",
      "selection_frequency_upper_95", "sign_stability",
      "sign_agreement_fraction", "stability_weight", "stable_estimate",
      "active_edge", "sign_flip_flag", "cv_rsq",
      "cv_predictive_above_null", "cv_preprocessing", "bootstrap_method",
      "n_bootstrap_requested", "n_bootstrap_success",
      "bootstrap_success_fraction", "min_bootstrap_success_fraction",
      "bootstrap_completion_adequate"
    ), colnames(condition_all))
    answer$stability_diagnostics <- condition_all[
      , stability_columns, drop = FALSE
    ]
    outdir <- call_args$outdir
    if (!is.null(outdir) && nzchar(as.character(outdir))) {
      .rc_mm_write_tsv_gz(
        answer$stability_diagnostics,
        file.path(outdir, "bootstrap_stability_diagnostics.tsv.gz")
      )
    }
  }
  answer$normalization_policy$cross_validation <- paste(
    "condition-stratified cell folds; condition centers and edge scales are",
    "estimated from each training fold and applied to its validation fold"
  )
  answer$normalization_policy$candidate_screen <- paste(
    "structural Pando candidates followed by condition-aware non-zero",
    "observability; outcome-based correlation threshold fixed at zero"
  )
  answer$normalization_policy$coefficient_parameterization <- paste(
    "elastic-net sparsity is applied directly to condition-specific theta;",
    "the global backbone is the cross-condition mean and deviations are derived"
  )
  answer
}
