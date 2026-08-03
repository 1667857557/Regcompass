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
  id_col <- if ("unit_id" %in% colnames(meta)) "unit_id" else "pool_id"
  unit_id <- as.character(meta[[id_col]])
  if (!setequal(unit_id, colnames(normalized))) {
    stop("Condition increment units and metadata differ.",
         call. = FALSE)
  }
  meta <- meta[match(colnames(normalized), unit_id), , drop = FALSE]
  row_meta <- rc_parse_microcompass_row_id(rownames(normalized))
  strata <- unique(meta[, c(condition_col, celltype_col), drop = FALSE])
  rows <- lapply(seq_len(nrow(strata)), function(i) {
    selected <- as.character(meta[[condition_col]]) ==
      as.character(strata[[condition_col]][[i]]) &
      as.character(meta[[celltype_col]]) ==
      as.character(strata[[celltype_col]][[i]])
    value <- matrixStats::rowMedians(
      normalized[, selected, drop = FALSE],
      na.rm = TRUE
    )
    value[is.nan(value)] <- NA_real_
    data.frame(
      row_id = rownames(normalized),
      reaction_id = row_meta$reaction_id,
      target_direction = row_meta$target_direction,
      medium_scenario = row_meta$medium_scenario,
      condition = as.character(strata[[condition_col]][[i]]),
      cell_type = as.character(strata[[celltype_col]][[i]]),
      median_condition_unique_increment_per_required_target_flux = value,
      n_metacells = sum(selected),
      inference_class = "metacell_statistical_unit_within_dataset",
      statistical_unit = "metacell",
      metacell_statistical_inference = TRUE,
      biological_replicate_inference = FALSE,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
