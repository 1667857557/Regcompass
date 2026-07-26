# Adapt meta-module provenance to the GRN evidence schema without changing GPR math.

.rc_summarize_supported_metabolic_genes_core <-
  .rc_summarize_supported_metabolic_genes

.rc_summarize_supported_metabolic_genes <- function(
    grn_result, metabolic_genes) {
  answer <- .rc_summarize_supported_metabolic_genes_core(
    grn_result = grn_result,
    metabolic_genes = metabolic_genes
  )
  multitask <- identical(
    grn_result$grn_mode %||% "legacy_condition_pando",
    "multitask_shared_backbone"
  )
  if (!multitask || !nrow(answer)) return(answer)

  significant <- grn_result$tf_peak_gene_significant
  key <- paste(significant$group_id, significant$target, sep = "\001")
  rows <- split(seq_len(nrow(significant)), key)
  diagnostics <- do.call(rbind, lapply(rows, function(index) {
    one <- significant[index, , drop = FALSE]
    finite_min <- function(x) {
      x <- suppressWarnings(as.numeric(x))
      x <- x[is.finite(x)]
      if (length(x)) min(x) else NA_real_
    }
    finite_max <- function(x) {
      x <- suppressWarnings(as.numeric(x))
      x <- x[is.finite(x)]
      if (length(x)) max(x) else NA_real_
    }
    data.frame(
      group_id = as.character(one$group_id[[1L]]),
      gene = toupper(as.character(one$target[[1L]])),
      min_selection_frequency = finite_min(one$selection_frequency),
      min_sign_stability = finite_min(one$sign_stability),
      max_abs_effective_estimate = finite_max(abs(one$effective_estimate)),
      any_sign_flip = any(one$sign_flip_flag %in% TRUE),
      stringsAsFactors = FALSE
    )
  }))
  answer <- merge(
    answer,
    diagnostics,
    by = c("group_id", "gene"),
    all.x = TRUE,
    sort = FALSE
  )
  answer$evidence_definition <-
    "stability_selected_multitask_tf_peak_gene_target"
  answer
}

.rc_build_condition_meta_modules_core <- .rc_build_condition_meta_modules
.rc_build_condition_meta_modules <- function(
    grn_result, gem, outdir, meta_module_args = list()) {
  answer <- .rc_build_condition_meta_modules_core(
    grn_result = grn_result,
    gem = gem,
    outdir = outdir,
    meta_module_args = meta_module_args
  )
  multitask <- identical(
    grn_result$grn_mode %||% "legacy_condition_pando",
    "multitask_shared_backbone"
  )
  if (multitask) {
    answer$core_definition <- paste(
      "complete GEM GPR branch contained in the stability-selected multitask",
      "target-gene set for one condition-by-cell-type group"
    )
    answer$analysis_group_unit <-
      "condition_x_celltype_stability_selected_multitask_metabolic_targets"
    answer$regulatory_gene_definition <- paste(
      "a metabolic target is condition-active when at least one TF-peak-target",
      "edge passes CV reliability, selection-frequency, sign-stability and",
      "effect-size thresholds; both positive and negative stable edges qualify"
    )
    if (is.data.frame(grn_result$condition_target_genes)) {
      answer$condition_target_genes <- grn_result$condition_target_genes
      .rc_mm_write_tsv_gz(
        answer$condition_target_genes,
        file.path(outdir, "condition_target_genes.tsv.gz")
      )
    }
    saveRDS(answer, file.path(outdir, "condition_meta_modules.rds"))
  }
  answer
}
