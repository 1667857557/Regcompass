.rc_condition_gene_regulatory_modifier_unfiltered <-
  .rc_condition_gene_regulatory_modifier

# Canonical correction: non-finite Pando R-squared is not regulatory evidence.
# Edges without finite R-squared are excluded; a target with no finite R-squared
# therefore receives a zero (untrusted) ATAC modifier rather than reliability 1.
.rc_condition_gene_regulatory_modifier_finite_rsq <- function(
    significant_edges, object, unit_meta,
    condition_col = "condition", celltype_col = "cell_type",
    atac_assay = "ATAC",
    target_genes = NULL,
    min_scale = 0.05) {
  edges <- significant_edges
  if (is.data.frame(edges) && nrow(edges)) {
    rsq <- if ("rsq" %in% colnames(edges)) {
      suppressWarnings(as.numeric(edges$rsq))
    } else {
      rep(NA_real_, nrow(edges))
    }
    edges <- edges[is.finite(rsq), , drop = FALSE]
  }
  out <- .rc_condition_gene_regulatory_modifier_unfiltered(
    significant_edges = edges,
    object = object,
    unit_meta = unit_meta,
    condition_col = condition_col,
    celltype_col = celltype_col,
    atac_assay = atac_assay,
    target_genes = target_genes,
    min_scale = min_scale
  )
  attr(out, "reliability_policy") <- paste(
    "only finite Pando R-squared values are trusted; targets without finite",
    "R-squared receive regulatory reliability zero"
  )
  out
}
.rc_condition_gene_regulatory_modifier <-
  .rc_condition_gene_regulatory_modifier_finite_rsq

.rc_feasibility_completion_metadata <- function(model_mode) {
  if (identical(model_mode, "meta_module_gem")) {
    return(list(
      feasibility_completion =
        "local_unconstrained_fastcore_then_global_union_medium_specific_fastcore",
      feasibility_completion_stages = list(
        local = paste(
          "condition x cell-type biological meta-modules are completed against",
          "an unconstrained shared FASTCC parent"
        ),
        global = paste(
          "the union model is rebuilt and add-only FASTCORE-completed separately",
          "for each shared medium scenario before scoring"
        )
      )
    ))
  }
  list(
    feasibility_completion = "not_applicable_full_gem",
    feasibility_completion_stages = list(
      local = "local meta-modules may be built upstream but are not used for full-GEM scoring",
      global = "the complete GEM is constrained by each shared medium without meta-module completion"
    )
  )
}

.rc_apply_corrected_result_metadata <- function(result) {
  completion <- .rc_feasibility_completion_metadata(result$model_mode)
  result$params$feasibility_completion <- completion$feasibility_completion
  result$params$feasibility_completion_stages <- completion$feasibility_completion_stages
  result$params$atac_tfidf_scope <- "cell_type_across_conditions"
  result$params$pando_normalization_policy <-
    "reuse Step 1 normalized RNA and cell-type-shared ATAC TF-IDF"
  result$params$zero_count_peak_policy <-
    "exclude before shared TF-IDF and within each Pando group"
  result$params$missing_expression_policy <-
    "unmeasured reaction expression is zero-filled and receives the same strict penalty as explicit zero"
  result$params$fragment_registration_policy <-
    "ignore and clear registered Signac fragments when fragment_files = FALSE"
  result
}

.rc_run_regcompass_uncorrected_metadata <- rc_run_regcompass
.rc_run_regcompass_v170 <- function(
    object, gem, outdir, pfm, genome,
    fragment_files = FALSE,
    sample_col = "sample_id",
    condition_col = "condition",
    celltype_col = "cell_type",
    rna_assay = "RNA",
    atac_assay = "ATAC",
    model_mode = c("meta_module_gem", "full_gem"),
    medium_scenarios = NULL,
    metacell_args = list(),
    layer1_args = list(),
    pando_args = list(),
    layer2_args = list(),
    upstream_workers = NULL,
    layer2_workers = NULL,
    parallel_backend = c("auto", "serial", "snow", "multicore"),
    species = c("auto", "human", "mouse")) {
  if (identical(fragment_files, FALSE) || is.null(fragment_files)) {
    object <- .rc_clear_signac_fragments(object, atac_assay = atac_assay)
  }
  .rc_condition_pool_design_summary(
    object@meta.data,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    strict_biological_defaults = FALSE
  )
  warning(
    paste(
      "Metacell-level scores are descriptive pseudo-observations and are not",
      "independent biological replicates."
    ),
    call. = FALSE
  )
  result <- .rc_run_regcompass_uncorrected_metadata(
    object = object, gem = gem, outdir = outdir, pfm = pfm, genome = genome,
    fragment_files = fragment_files,
    sample_col = sample_col, condition_col = condition_col,
    celltype_col = celltype_col, rna_assay = rna_assay,
    atac_assay = atac_assay, model_mode = model_mode,
    medium_scenarios = medium_scenarios,
    metacell_args = metacell_args, layer1_args = layer1_args,
    pando_args = pando_args, layer2_args = layer2_args,
    upstream_workers = upstream_workers, layer2_workers = layer2_workers,
    parallel_backend = parallel_backend, species = species
  )
  result <- .rc_apply_corrected_result_metadata(result)
  saveRDS(result, file.path(outdir, "regcompass_result.rds"))
  rc_export_microcompass(result$microcompass, outdir)
  result
}
rc_run_regcompass <- .rc_run_regcompass_v170

.rc_regcompass_step_results_uncorrected_metadata <- rc_regcompass_step_results
.rc_regcompass_step_results_v170 <- function(
    metacells, meta_modules, layer1, layer2, gem, outdir,
    species = c("auto", "human", "mouse")) {
  result <- .rc_regcompass_step_results_uncorrected_metadata(
    metacells = metacells,
    meta_modules = meta_modules,
    layer1 = layer1,
    layer2 = layer2,
    gem = gem,
    outdir = outdir,
    species = species
  )
  result <- .rc_apply_corrected_result_metadata(result)
  saveRDS(result, file.path(outdir, "regcompass_result.rds"))
  result
}
rc_regcompass_step_results <- .rc_regcompass_step_results_v170
