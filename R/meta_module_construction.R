.rc_build_condition_meta_modules <- function(
    grn_result, gem, outdir, layer1_args = list()) {
  if (!is.list(grn_result) ||
      !is.data.frame(grn_result$tf_peak_gene_significant)) {
    stop("`grn_result` is not a valid single-cell GRN result.", call. = FALSE)
  }
  obsolete <- intersect(
    names(layer1_args), c("local_fastcore", "local_fastcore_args")
  )
  if (length(obsolete)) {
    stop(
      paste0(
        "Local FASTCORE was removed. Delete obsolete `layer1_args` fields: ",
        paste(obsolete, collapse = ", "),
        ". Configure the single medium-specific global FASTCORE through ",
        "`layer2_args$model_params`."
      ),
      call. = FALSE
    )
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  group_cols <- grn_result$group_cols
  display_cols <- c("group_id", group_cols)
  module_cols <- unique(c(display_cols, "sample_id", "module_id"))
  metabolic_genes <- gem$metabolic_genes %||%
    rc_metabolic_gpr_genes(gem$gpr_table)
  sig <- grn_result$tf_peak_gene_significant
  sig$sample_id <- sig$group_id
  projection <- rc_project_metabolic_grn(
    sig,
    metabolic_genes = metabolic_genes,
    top_k = layer1_args$top_k_neighbors %||% 5L,
    min_shared_tfs = layer1_args$min_shared_tfs %||% 1L,
    min_tf_jaccard = layer1_args$min_tf_jaccard %||% 0,
    max_targets_per_tf = layer1_args$max_targets_per_tf %||% 200L,
    include_direct_metabolic_tf = TRUE
  )
  group_meta <- unique(
    grn_result$sample_status[, display_cols, drop = FALSE]
  )
  group_meta$analysis_unit_id <- group_meta$group_id
  projection$nodes <- .rc_remap_projection_metadata(
    projection$nodes, group_meta, "analysis_unit_id", display_cols
  )
  projection$edges <- .rc_remap_projection_metadata(
    projection$edges, group_meta, "analysis_unit_id", display_cols
  )
  core <- rc_map_meta_module_core_reactions(
    projection$nodes, gem$gpr_table
  )
  if (nrow(core)) {
    core <- merge(
      core,
      unique(projection$nodes[, module_cols, drop = FALSE]),
      by = c("sample_id", "module_id"),
      all.x = TRUE,
      sort = FALSE
    )
    core <- core[, c(
      display_cols, setdiff(colnames(core), display_cols)
    ), drop = FALSE]
  }
  expanded <- rc_expand_meta_module_reactions(
    gem,
    core,
    subsystem_table = layer1_args$subsystem_table %||% NULL,
    expansion_mode = layer1_args$expansion_mode %||% "ordered_once"
  )
  if (nrow(expanded$reaction_membership)) {
    expanded$reaction_membership <- merge(
      expanded$reaction_membership,
      unique(core[, module_cols, drop = FALSE]),
      by = c("sample_id", "module_id"),
      all.x = TRUE,
      sort = FALSE
    )
    expanded$reaction_membership <- expanded$reaction_membership[, c(
      display_cols,
      setdiff(colnames(expanded$reaction_membership), display_cols)
    ), drop = FALSE]
  }
  out <- c(grn_result, list(
    metabolic_gene_nodes = projection$nodes,
    metabolic_gene_edges = projection$edges,
    core_gene_reaction = core,
    reaction_membership = expanded$reaction_membership,
    biological_reaction_membership = expanded$reaction_membership,
    meta_module_summary = expanded$summary,
    crossref_maps = expanded$crossref_maps,
    analysis_group_unit = "condition_x_celltype_single_cell_grn",
    feasibility_completion = "none_at_meta_module_stage"
  ))
  .rc_mm_write_tsv_gz(
    projection$nodes, file.path(outdir, "metabolic_gene_nodes.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    projection$edges, file.path(outdir, "metabolic_gene_edges.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    core, file.path(outdir, "core_gene_reaction.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    expanded$reaction_membership,
    file.path(outdir, "meta_module_reactions.tsv.gz")
  )
  saveRDS(out, file.path(outdir, "condition_meta_modules.rds"))
  out
}
