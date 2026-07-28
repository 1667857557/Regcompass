# Clarify complete-branch and reaction-core semantics after the core mapper.

.rc_map_meta_module_core_reactions_raw <- rc_map_meta_module_core_reactions

rc_map_meta_module_core_reactions <- function(gene_nodes, gpr_table) {
  out <- .rc_map_meta_module_core_reactions_raw(gene_nodes, gpr_table)
  if (!is.data.frame(out) || !nrow(out)) {
    if (is.data.frame(out) && !"reaction_is_core" %in% colnames(out)) {
      out$reaction_is_core <- logical(nrow(out))
    }
    return(out)
  }
  required <- c(
    "group_id", "module_id", "reaction_id", "and_group_id",
    "group_complete"
  )
  missing <- setdiff(required, colnames(out))
  if (length(missing)) {
    stop(
      "GPR core mapping output is missing columns: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  reaction_key <- paste(
    out$group_id, out$module_id, out$reaction_id, sep = "\001"
  )
  reaction_complete <- tapply(
    out$group_complete %in% TRUE,
    reaction_key,
    any
  )
  out$reaction_is_core <- unname(reaction_complete[reaction_key])
  out$is_core <- out$group_complete %in% TRUE
  out$is_partial_candidate <- !out$group_complete
  out$inclusion_stage <- ifelse(
    out$group_complete,
    "core_complete_gpr",
    ifelse(
      out$reaction_is_core,
      "alternative_incomplete_gpr_branch",
      "partial_gpr_candidate"
    )
  )
  out
}
