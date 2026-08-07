# Internal runtime RDS storage policy.
#
# Package functions call `saveRDS()` without namespace qualification. Defining
# the function in the package namespace centralizes checkpoint storage rules:
#
# 1. runtime RDS files use explicit gzip compression while retaining the
#    conventional `.rds` extension;
# 2. one authoritative checkpoint is retained per workflow stage;
# 3. large RDS files already embedded in a stage checkpoint are not written
#    again as side artifacts;
# 4. expensive structural GEM caches remain external because stage checkpoints
#    contain only their references and checksums;
# 5. stage checkpoints themselves are compacted only where they contain exact
#    aliases or complete upstream-stage payloads that are not needed to reload
#    the stage for downstream RegCompass steps.

.rc_duplicate_runtime_rds <- function(file) {
  if (!is.character(file) || length(file) != 1L ||
      is.na(file) || !nzchar(file)) {
    return(FALSE)
  }

  name <- basename(file)
  duplicate_names <- c(
    # Stage 1: complete Pando/GRN objects already live in step_grn.rds.
    "single_cell_grn.rds",
    "pando_condition_grn_fits.rds",
    "condition_grn_fit.rds",

    # Stage 2: the canonical metacell object lives once in step_metacells.rds.
    "rna_counts.rds",
    "atac_counts.rds",
    "metacell_object.rds",
    "merged_metacell_object.rds",
    "condition_metacell_cache_contract.rds",

    # Stage 3: merged modules are embedded in step_meta_modules.rds. The
    # condition-level module file is intentionally NOT listed because Stage 3
    # stores it externally and keeps only a checksum-validated reference.
    "merged_meta_modules.rds",

    # Stage 5: these matrices/contracts are already fields of step_layer2.rds.
    "strict_score_matrix.rds",
    "strict_penalty_matrix.rds",
    "vmax_matrix.rds",
    "feasible_matrix.rds",
    "evaluated_matrix.rds",
    "penalty_components.rds",
    "run_parameters.rds",
    "model_completion_contract.rds",

    # Stage 6: comparison summaries are already fields of regcompass_result.rds.
    "step_comparison.rds",

    # One-shot convenience files duplicate stage-owned information.
    "00_model_info.rds",
    "00_medium_scenarios.rds"
  )

  if (name %in% duplicate_names) {
    return(TRUE)
  }

  # Standard Pando previously wrote one complete GRN object per cell type.
  if (grepl("^standard_pando_.*\\.rds$", name)) {
    return(TRUE)
  }

  # rc_run_regcompass() used to write the final result a second time at the
  # workflow root after Stage 6 had already written the canonical copy.
  if (identical(name, "regcompass_result.rds")) {
    canonical <- file.path(
      dirname(file), "06_results", "regcompass_result.rds"
    )
    if (file.exists(canonical)) {
      return(TRUE)
    }
  }

  FALSE
}

.rc_compact_step_grn_storage <- function(object) {
  if (!is.list(object) || !is.list(object$grn_result)) {
    return(object)
  }

  grn <- object$grn_result

  # In the parallel condition-GRN contract the named object map is canonical.
  # For a single condition-GRN cell type `pando_grn_data` is an exact second
  # reference to the same GRNData object and is unnecessary after reload.
  object_map <- grn$pando_grn_data_by_cell_type
  if (is.list(object_map) && length(object_map)) {
    grn$pando_grn_data <- NULL
  }

  # These historical aliases are exact copies of the canonical condition-edge
  # tables. Downstream RegCompass code uses tf_peak_gene_condition[_all].
  grn$tf_peak_gene_condition_effect_all <- NULL
  grn$tf_peak_gene_condition_effect <- NULL

  object$grn_result <- grn
  object
}

.rc_compact_step_metacells_storage <- function(object) {
  if (!is.list(object) || !is.list(object$pooled)) {
    return(object)
  }

  pooled <- object$pooled

  # `object$metacell_object` is the single canonical Seurat metacell object.
  # The pooled copies and explicit RNA/ATAC count matrices are already contained
  # in that object's assays and can otherwise multiply Stage 2 disk usage.
  for (field in c(
    "metacell_objects", "metacell_object", "rna_counts", "atac_counts"
  )) {
    pooled[[field]] <- NULL
  }

  object$pooled <- pooled
  object
}

.rc_compact_step_layer2_storage <- function(object) {
  if (!is.list(object)) {
    return(object)
  }

  # The RNA-only route previously retained another near-complete microCOMPASS
  # result even though its required matrices are already stored at the top level
  # as penalty_rna_only and score_rna_only and share the primary model cache.
  object$comparison_paths <- NULL
  object
}

.rc_compact_final_result_storage <- function(object) {
  if (!is.list(object)) {
    return(object)
  }

  upstream_fields <- intersect(c(
    "grn",
    "metacells",
    "layer1",
    "condition_grn_meta_modules",
    "merged_grn_meta_modules",
    "grn_meta_modules",
    "microcompass"
  ), names(object))

  for (field in upstream_fields) {
    object[[field]] <- NULL
  }

  existing <- object$storage_contract
  if (!is.list(existing)) {
    existing <- list()
  }
  object$storage_contract <- utils::modifyList(existing, list(
    schema_version = "regcompass_compact_final_result_storage_v1",
    upstream_stage_payloads_embedded = FALSE,
    omitted_upstream_fields = upstream_fields,
    final_analysis_tables_retained = intersect(c(
      "reaction_catalog",
      "reaction_evidence",
      "reaction_comparison_by_metacell",
      "reaction_ranking",
      "condition_summary",
      "condition_contrast",
      "rna_only_control_summary",
      "rna_only_control_contrast"
    ), names(object)),
    upstream_policy = paste(
      "load the corresponding step_grn/step_metacells/step_meta_modules/",
      "step_layer1/step_layer2 checkpoint when upstream internals are needed"
    )
  ))

  object
}

.rc_compact_runtime_rds <- function(object, file) {
  if (!is.character(file) || length(file) != 1L ||
      is.na(file) || !nzchar(file)) {
    return(object)
  }

  switch(
    basename(file),
    step_grn.rds = .rc_compact_step_grn_storage(object),
    step_metacells.rds = .rc_compact_step_metacells_storage(object),
    step_layer2.rds = .rc_compact_step_layer2_storage(object),
    regcompass_result.rds = .rc_compact_final_result_storage(object),
    object
  )
}

saveRDS <- function(
    object, file = "", ascii = FALSE, version = NULL,
    compress = TRUE, refhook = NULL) {
  is_path <- is.character(file) && length(file) == 1L &&
    !is.na(file) && nzchar(file)

  if (is_path && .rc_duplicate_runtime_rds(file)) {
    return(invisible(NULL))
  }

  if (is_path) {
    object <- .rc_compact_runtime_rds(object, file)
  }

  resolved_compression <- if (is_path) "gzip" else compress
  base::saveRDS(
    object = object,
    file = file,
    ascii = ascii,
    version = version,
    compress = resolved_compression,
    refhook = refhook
  )
}
