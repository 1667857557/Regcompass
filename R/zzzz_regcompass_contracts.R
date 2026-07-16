.rc_run_regcompass_stratum_balanced_v15 <- .rc_run_regcompass_stratum
rc_run_regcompass_balanced_v15 <- rc_run_regcompass

.rc_run_regcompass_stratum <- function(object, group_id, group_cols, gem, outdir,
                                        pfm, genome, fragment_files = NULL,
                                        sample_col = "sample_id",
                                        condition_col = "condition",
                                        celltype_col = "cell_type",
                                        rna_assay = "RNA",
                                        atac_assay = "ATAC",
                                        metacell_args = list(),
                                        layer1_args = list(),
                                        pando_args = list()) {
  result <- .rc_run_regcompass_stratum_balanced_v15(
    object = object,
    group_id = group_id,
    group_cols = group_cols,
    gem = gem,
    outdir = outdir,
    pfm = pfm,
    genome = genome,
    fragment_files = fragment_files,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    metacell_args = metacell_args,
    layer1_args = layer1_args,
    pando_args = pando_args
  )
  if (identical(result$status, "ok")) {
    artifact <- readRDS(result$artifact_file)
    artifact$schema_version <- "regcompass_stratum_v2"
    artifact$architecture_version <-
      "local_fastcore_sample_balanced_global_calibration_v1"
    saveRDS(artifact, result$artifact_file)
  }
  result
}

rc_run_regcompass <- function(object, gem, outdir, pfm, genome,
                               fragment_files = NULL,
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
                               parallel_backend = c("auto", "serial", "snow", "multicore")) {
  result <- rc_run_regcompass_balanced_v15(
    object = object,
    gem = gem,
    outdir = outdir,
    pfm = pfm,
    genome = genome,
    fragment_files = fragment_files,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    model_mode = model_mode,
    medium_scenarios = medium_scenarios,
    metacell_args = metacell_args,
    layer1_args = layer1_args,
    pando_args = pando_args,
    layer2_args = layer2_args,
    upstream_workers = upstream_workers,
    layer2_workers = layer2_workers,
    parallel_backend = parallel_backend
  )
  result$params$capacity_calibration_scope <-
    result$layer1$capacity_calibration_scope
  result$params$sample_balanced_gene_score <-
    isTRUE(result$layer1$calibration_params$sample_balance)
  result$params$sample_balanced_q95 <-
    isTRUE(result$layer1$calibration_params$sample_balance)
  result$params$expression_batch_correction <-
    result$layer1$calibration_params$expression_batch_correction
  result$params$local_fastcore_before_global_union <-
    identical(
      result$grn_meta_modules$global_union_source,
      "deduplicated_local_fastcore_completed_meta_modules"
    )
  result$params$global_fastcore_repair <-
    "conditional_on_global_core_directional_incompleteness"
  saveRDS(result, file.path(outdir, "regcompass_global_metacell_result.rds"))
  result
}
