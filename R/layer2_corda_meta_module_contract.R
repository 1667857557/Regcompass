# Resolve Stage 3 biological reaction membership for CORDA2 evidence mapping.

.rc_meta_module_reaction_membership <- function(meta_modules) {
  if (!is.list(meta_modules)) {
    stop("`meta_modules` must be a Stage 3 result or merged catalogue.",
         call. = FALSE)
  }

  roots <- list(meta_modules$merged_modules, meta_modules)
  roots <- roots[vapply(roots, is.list, logical(1))]
  fields <- c(
    "merged_reaction_membership",
    "reaction_membership",
    "biological_reaction_membership"
  )
  celltype_candidates <- unique(c(
    "cell_type",
    as.character(meta_modules$workflow_params$celltype_col),
    as.character(meta_modules$celltype_col),
    as.character(meta_modules$merged_modules$celltype_col)
  ))
  celltype_candidates <- celltype_candidates[
    !is.na(celltype_candidates) & nzchar(celltype_candidates)
  ]

  for (root in roots) {
    root_celltype <- as.character(root$celltype_col)
    root_candidates <- unique(c(celltype_candidates, root_celltype))
    root_candidates <- root_candidates[
      !is.na(root_candidates) & nzchar(root_candidates)
    ]
    for (field in fields) {
      tab <- root[[field]]
      if (!is.data.frame(tab) || !nrow(tab) ||
          !"reaction_id" %in% colnames(tab)) {
        next
      }
      celltype_col <- intersect(root_candidates, colnames(tab))
      if (!length(celltype_col)) next
      celltype_col <- celltype_col[[1L]]

      tab$reaction_id <- trimws(as.character(tab$reaction_id))
      tab$cell_type <- trimws(as.character(tab[[celltype_col]]))
      keep <- !is.na(tab$reaction_id) & nzchar(tab$reaction_id) &
        !is.na(tab$cell_type) & nzchar(tab$cell_type)
      tab <- tab[keep, , drop = FALSE]
      if (!nrow(tab)) next

      key <- paste(tab$cell_type, tab$reaction_id, sep = "\001")
      tab <- tab[!duplicated(key), , drop = FALSE]
      rownames(tab) <- NULL
      return(tab)
    }
  }

  stop(
    "Stage 3 does not contain non-empty merged reaction membership with ",
    "reaction_id and cell-type columns required for CORDA2 evidence mapping.",
    call. = FALSE
  )
}
