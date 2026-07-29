# Merge biological meta-modules without constructing a GEM.

.rc_merge_meta_module_catalogue <- function(condition_modules) {
  names_to_merge <- c(
    "condition_fit_status", "tf_peak_gene_condition_all",
    "tf_peak_gene_condition",
    "supported_metabolic_genes", "core_gene_reaction",
    "reaction_membership", "meta_module_summary"
  )
  out <- lapply(names_to_merge, function(name) {
    value <- condition_modules[[name]]
    if (is.data.frame(value)) value else data.frame()
  })
  names(out) <- names_to_merge

  core <- out$core_gene_reaction
  if ("is_core" %in% colnames(core)) {
    core <- core[core$is_core %in% TRUE, , drop = FALSE]
  }
  core_ids <- unique(as.character(core$reaction_id))
  core_ids <- core_ids[!is.na(core_ids) & nzchar(core_ids)]

  biological <- out$reaction_membership
  biological_ids <- unique(as.character(biological$reaction_id))
  biological_ids <- biological_ids[
    !is.na(biological_ids) & nzchar(biological_ids)
  ]
  if (!length(core_ids) || !length(biological_ids)) {
    stop(
      "No merged biological meta-module reactions were produced.",
      call. = FALSE
    )
  }
  if (length(setdiff(core_ids, biological_ids))) {
    stop(
      "Merged core reactions are missing from biological membership.",
      call. = FALSE
    )
  }

  out$biological_reaction_membership <- biological
  out$merged_core_reactions <- data.frame(
    catalogue_id = "MERGED_META_MODULES",
    reaction_id = core_ids,
    is_core = TRUE,
    stringsAsFactors = FALSE
  )
  out$merged_reaction_membership <- data.frame(
    catalogue_id = "MERGED_META_MODULES",
    reaction_id = biological_ids,
    is_core = biological_ids %in% core_ids,
    inclusion_stage = ifelse(
      biological_ids %in% core_ids,
      "merged_meta_module_core",
      "merged_meta_module_biological_member"
    ),
    stringsAsFactors = FALSE
  )
  out$schema_version <- "regcompass_merged_meta_modules_v2"
  out$source_group_ids <- if (
    "group_id" %in% colnames(out$condition_fit_status)
  ) {
    unique(as.character(out$condition_fit_status$group_id))
  } else {
    character()
  }
  out$core_definition <-
    "condition_celltype_active_pando_targets_complete_gpr"
  out$merge_source <- "deduplicated_biological_meta_module_reactions"
  out$is_gem <- FALSE
  out$fastcore_applied <- FALSE
  out
}
