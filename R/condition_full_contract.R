# Common-dictionary condition result helpers.

# The Pando fit and extraction validators are defined by condition_grn_contract.R.
# This file intentionally contains no OOF, common-support, structural-zero or
# projection-mask overrides.

.rc_condition_increment_summary <- function(
    microcompass, condition_col, celltype_col) {
  increment <- as.matrix(
    microcompass$penalty_condition_unique_increment
  )
  vmax <- as.matrix(microcompass$vmax)
  if (!identical(dimnames(increment), dimnames(vmax))) {
    stop("Condition-unique increment and vmax matrices are not aligned.",
         call. = FALSE)
  }
  omega <- microcompass$params$omega %||% 0.95
  normalized <- increment / (omega * vmax)
  normalized[
    !is.finite(increment) | !is.finite(vmax) | vmax <= 0
  ] <- NA_real_
  meta <- microcompass$unit_meta
  if (!is.data.frame(meta) ||
      !all(c(condition_col, celltype_col) %in% colnames(meta))) {
    stop("Condition increment metadata lack condition/cell-type columns.",
         call. = FALSE)
  }
  id_col <- if ("unit_id" %in% colnames(meta)) {
    "unit_id"
  } else if ("pool_id" %in% colnames(meta)) {
    "pool_id"
  } else {
    stop("Condition increment metadata lack unit_id/pool_id.", call. = FALSE)
  }
  unit_id <- trimws(as.character(meta[[id_col]]))
  if (anyNA(unit_id) || any(!nzchar(unit_id)) || anyDuplicated(unit_id) ||
      !setequal(unit_id, colnames(normalized))) {
    stop("Condition increment units and metadata differ.",
         call. = FALSE)
  }
  meta <- meta[match(colnames(normalized), unit_id), , drop = FALSE]
  row_meta <- rc_parse_microcompass_row_id(rownames(normalized))
  row_meta$row_id <- rownames(normalized)
  strata <- unique(meta[, c(condition_col, celltype_col), drop = FALSE])
  rows <- lapply(seq_len(nrow(strata)), function(i) {
    condition <- as.character(strata[[condition_col]][[i]])
    cell_type <- as.character(strata[[celltype_col]][[i]])
    selected <- as.character(meta[[condition_col]]) == condition &
      as.character(meta[[celltype_col]]) == cell_type
    row_keep <- is.na(row_meta$cell_type) | row_meta$cell_type == cell_type
    if (!any(row_keep)) return(NULL)
    value <- matrixStats::rowMedians(
      normalized[row_keep, selected, drop = FALSE],
      na.rm = TRUE
    )
    value[is.nan(value)] <- NA_real_
    data.frame(
      row_id = row_meta$row_id[row_keep],
      reaction_id = row_meta$reaction_id[row_keep],
      target_direction = row_meta$target_direction[row_keep],
      medium_scenario = row_meta$medium_scenario[row_keep],
      condition = condition,
      cell_type = cell_type,
      median_condition_unique_increment_per_required_target_flux = value,
      n_metacells = sum(selected),
      inference_class = "metacell_statistical_unit_within_dataset",
      statistical_unit = "metacell",
      metacell_statistical_inference = TRUE,
      biological_replicate_inference = FALSE,
      structural_scope = "same_celltype_conditions_on_celltype_medium_union_gem",
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(data.frame())
  answer <- do.call(rbind, rows)
  rownames(answer) <- NULL
  answer
}
