# Integrated RegCompass v1.3 runner.

#' Run the integrated RegCompass workflow with one of two GEM modes
#' @export
rc_run_regcompass <- function(object, gem, outdir, pfm, genome,
                               fragment_files = NULL,
                               sample_col = "sample_id",
                               condition_col = "condition",
                               celltype_col = "cell_type",
                               rna_assay = "RNA",
                               atac_assay = "ATAC",
                               model_mode = c(
                                 "meta_module_gem", "full_gem"
                               ),
                               medium_scenarios = NULL,
                               metacell_args = list(),
                               pando_args = list(),
                               layer2_args = list()) {
  model_mode <- match.arg(model_mode)
  if (is.null(gem$gpr_table)) {
    stop("`gem` must contain `gpr_table`.", call. = FALSE)
  }
  reserved_layer2 <- intersect(
    names(layer2_args),
    c(
      "layer1", "gem", "mode", "reaction_membership",
      "core_reactions", "medium_scenarios", "sample_col",
      "condition_col", "celltype_col"
    )
  )
  if (length(reserved_layer2)) {
    stop(
      paste0(
        "`layer2_args` cannot override integrated workflow fields: ",
        paste(reserved_layer2, collapse = ", "), "."
      ),
      call. = FALSE
    )
  }

  layer1_defaults <- list(
    object = object,
    gpr_table = gem$gpr_table,
    outdir = outdir,
    fragment_files = fragment_files,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay
  )
  layer1_defaults[names(metacell_args)] <- NULL
  result <- do.call(
    rc_run_regcompass_multiome_metacell,
    c(layer1_defaults, metacell_args)
  )

  retained <- as.character(result$metacell_meta$metacell_id)
  metacell_object <- rc_load_metacell_object_from_run(
    outdir,
    retained_metacell_ids = retained,
    rna_assay = rna_assay,
    atac_assay = atac_assay
  )
  single_cell_genes <- rownames(
    .rc_get_assay_counts(object, rna_assay)
  )
  pando_defaults <- list(
    metacell_object = metacell_object,
    gem = gem,
    outdir = file.path(outdir, "04_pando_meta_modules"),
    pfm = pfm,
    genome = genome,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    single_cell_genes = single_cell_genes,
    rna_assay = rna_assay,
    atac_assay = atac_assay
  )
  pando_defaults[names(pando_args)] <- NULL
  result$grn_meta_modules <- do.call(
    rc_run_pando_meta_modules,
    c(pando_defaults, pando_args)
  )

  core_gene_reaction <- result$grn_meta_modules$core_gene_reaction
  if ("is_core" %in% colnames(core_gene_reaction)) {
    core_gene_reaction <- core_gene_reaction[
      core_gene_reaction$is_core %in% TRUE,
      , drop = FALSE
    ]
  }
  default_targets <- unique(core_gene_reaction$reaction_id)
  layer2_defaults <- list(
    layer1 = result,
    gem = gem,
    target_reactions = default_targets,
    medium_scenarios = medium_scenarios,
    mode = model_mode,
    reaction_membership = if (
      identical(model_mode, "meta_module_gem")
    ) {
      result$grn_meta_modules$reaction_membership
    } else {
      NULL
    },
    core_reactions = if (
      identical(model_mode, "meta_module_gem")
    ) {
      result$grn_meta_modules$core_gene_reaction
    } else {
      NULL
    },
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col
  )
  layer2_defaults[names(layer2_args)] <- NULL
  result$microcompass <- do.call(
    rc_run_microcompass,
    c(layer2_defaults, layer2_args)
  )
  result$schema_version <- "regcompass_v1.3"
  result$model_mode <- model_mode
  saveRDS(
    result,
    file.path(outdir, "regcompass_v1.3_result.rds")
  )
  result
}
