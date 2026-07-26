# Compact final-result helpers. Detailed intermediate objects remain in their
# stage checkpoints instead of being duplicated inside regcompass_result.rds.

.rc_result_select_columns <- function(x, columns) {
  if (!is.data.frame(x)) return(data.frame())
  x[, intersect(columns, colnames(x)), drop = FALSE]
}

.rc_result_unique <- function(x, columns = colnames(x)) {
  x <- .rc_result_select_columns(x, columns)
  if (!nrow(x)) return(x)
  rownames(x) <- NULL
  unique(x)
}

.rc_compact_active_edges <- function(grn_result, condition_col, celltype_col) {
  edges <- grn_result$tf_peak_gene_significant
  columns <- c(
    "edge_id", condition_col, celltype_col, "tf", "atac_feature_id",
    "region", "target", "global_estimate", "condition_deviation",
    "effective_estimate", "selection_frequency", "sign_stability",
    "stable_estimate", "sign_flip_flag", "cv_rsq",
    "bootstrap_success_fraction", "active_edge", "evidence_type"
  )
  .rc_result_unique(edges, columns)
}

.rc_compact_condition_targets <- function(grn_result, condition_col, celltype_col) {
  targets <- grn_result$condition_target_genes
  columns <- c(
    "group_id", condition_col, celltype_col, "target", "n_active_edges",
    "n_positive_edges", "n_negative_edges", "max_selection_frequency",
    "max_sign_stability", "max_abs_effective_estimate",
    "max_abs_stable_estimate", "best_cv_rsq"
  )
  .rc_result_unique(targets, columns)
}

.rc_compact_core_reactions <- function(condition_modules) {
  core <- condition_modules$core_gene_reaction
  if (!is.data.frame(core) || !nrow(core)) return(data.frame())
  if ("reaction_is_core" %in% colnames(core)) {
    core <- core[core$reaction_is_core %in% TRUE, , drop = FALSE]
  } else if ("is_core" %in% colnames(core)) {
    core <- core[core$is_core %in% TRUE, , drop = FALSE]
  }
  columns <- c(
    "group_id", "module_id", "condition", "cell_type", "reaction_id",
    "reaction_is_core", "n_complete_gpr_branches", "inclusion_stage"
  )
  .rc_result_unique(core, columns)
}

.rc_compact_meta_module_summary <- function(condition_modules) {
  .rc_result_unique(
    condition_modules$meta_module_summary,
    c(
      "group_id", "module_id", "condition", "cell_type",
      "n_supported_genes", "n_core_reactions", "n_module_reactions"
    )
  )
}

.rc_compact_reaction_catalog <- function(catalog) {
  .rc_result_unique(catalog, c(
    "reaction_id", "reaction_name", "subsystem", "reaction_role",
    "lower_bound", "upper_bound", "reversible", "model_formula",
    "genes", "gpr_rule", "n_gpr_genes", "has_gpr",
    "kegg_reaction_id", "reactome_reaction_id", "rhea_reaction_id",
    "rhea_master_id"
  ))
}

.rc_compact_reaction_evidence <- function(evidence) {
  .rc_result_unique(evidence, c(
    "reaction_id", "condition", "cell_type", "evidence_class",
    "rna_capacity_median", "multiome_capacity_median",
    "delta_capacity_median", "n_metacells", "evidence_available"
  ))
}

.rc_result_table_manifest <- function(tables) {
  rows <- lapply(names(tables), function(name) {
    value <- tables[[name]]
    data.frame(
      table = name,
      n_rows = if (is.data.frame(value)) nrow(value) else NA_integer_,
      n_columns = if (is.data.frame(value)) ncol(value) else NA_integer_,
      role = switch(
        name,
        reaction_ranking = "primary within-condition reaction ranking",
        condition_contrast = "primary between-condition descriptive contrast",
        active_regulatory_edges = "bootstrap-active TF-peak-target evidence",
        condition_target_genes = "condition-specific regulated metabolic genes",
        core_reactions = "condition-specific complete-GPR reaction cores",
        meta_module_summary = "module-level counts only",
        grn_metacell_group_coverage = "GRN/metacell coverage audit",
        reaction_catalog = "unique annotations for scored reactions",
        reaction_evidence = "compact condition-by-cell-type evidence class",
        "compact analysis table"
      ),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.rc_result_intermediate_policy <- function() {
  data.frame(
    object = c(
      "grn candidate/all-coefficient tables",
      "metacell count matrices and Seurat object",
      "Layer 1 gene/reaction matrices",
      "full condition and merged meta-module membership",
      "Layer 2 auxiliary diagnostics"
    ),
    final_result = "not embedded",
    source = c(
      "Stage 1 step_grn.rds and TSV files",
      "Stage 2 step_metacells.rds and merged_metacell_object.rds",
      "Stage 4 step_layer1.rds",
      "Stage 3 step_meta_modules.rds and TSV files",
      "Stage 5 step_layer2.rds"
    ),
    stringsAsFactors = FALSE
  )
}
