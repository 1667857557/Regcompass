# Shared helpers for reaction annotation and evidence tables.

.rc_ra_nonempty <- function(x) {
  x <- as.character(x)
  !is.na(x) & nzchar(trimws(x))
}

.rc_ra_first_col <- function(x, candidates) {
  hit <- candidates[candidates %in% colnames(x)]
  if (length(hit)) hit[[1L]] else NULL
}

.rc_ra_collapse <- function(x, sep = ";") {
  x <- trimws(as.character(x))
  x <- unique(x[.rc_ra_nonempty(x)])
  if (length(x)) paste(x, collapse = sep) else NA_character_
}

.rc_ra_split <- function(x) {
  if (length(x) != 1L || is.na(x) || !nzchar(trimws(x))) return(character())
  out <- trimws(unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE))
  unique(out[nzchar(out)])
}

.rc_ra_infer_group_columns <- function(meta, condition_col, celltype_col) {
  if (!is.data.frame(meta)) {
    stop("Layer 1 `unit_meta` must be a data frame.", call. = FALSE)
  }
  if (is.null(condition_col)) {
    condition_col <- .rc_ra_first_col(
      meta, c("condition", "dataset", "Group", "group", "treatment")
    )
  }
  if (is.null(celltype_col)) {
    celltype_col <- .rc_ra_first_col(
      meta, c("cell_type", "celltype", "epithelial_or_stem", "CellType")
    )
  }
  if (is.null(condition_col) || is.null(celltype_col) ||
      !condition_col %in% colnames(meta) || !celltype_col %in% colnames(meta)) {
    stop(
      "Could not identify condition and cell-type columns in Layer 1 metadata.",
      call. = FALSE
    )
  }
  list(condition_col = condition_col, celltype_col = celltype_col)
}

.rc_ra_evidence_index <- function(evidence) {
  if (!is.data.frame(evidence) || !nrow(evidence)) return(character())
  paste(
    evidence$reaction_id,
    evidence$condition,
    evidence$cell_type,
    sep = "\001"
  )
}
