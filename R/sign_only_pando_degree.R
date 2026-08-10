# Target-degree normalization for the canonical sign-only Pando projection.

.rc_sign_only_target_degree <- function(
    grn_result, unit_meta, genes, condition_col, celltype_col) {
  genes <- as.character(genes)
  units <- as.character(unit_meta$unit_id)
  degree <- matrix(
    NA_real_, length(genes), length(units),
    dimnames = list(genes, units)
  )
  edge <- grn_result$tf_peak_gene_condition
  required <- c("target", condition_col, celltype_col)
  if (!is.data.frame(edge) || !nrow(edge)) return(degree)
  if (!all(required %in% colnames(edge))) {
    stop(
      "Sign-only target-degree normalization requires target, condition and ",
      "cell-type columns on active Pando edges.", call. = FALSE
    )
  }
  if (!all(c("unit_id", condition_col, celltype_col) %in% colnames(unit_meta))) {
    stop(
      "Sign-only target-degree normalization requires aligned metacell metadata.",
      call. = FALSE
    )
  }

  edge_target <- tolower(trimws(as.character(edge$target)))
  edge_condition <- as.character(edge[[condition_col]])
  edge_celltype <- as.character(edge[[celltype_col]])
  unit_condition <- as.character(unit_meta[[condition_col]])
  unit_celltype <- as.character(unit_meta[[celltype_col]])

  strata <- unique(data.frame(
    condition = unit_condition,
    cell_type = unit_celltype,
    stringsAsFactors = FALSE
  ))
  for (i in seq_len(nrow(strata))) {
    condition <- strata$condition[[i]]
    cell_type <- strata$cell_type[[i]]
    selected_units <- units[
      unit_condition == condition & unit_celltype == cell_type
    ]
    selected_edges <-
      edge_condition == condition & edge_celltype == cell_type
    if (!length(selected_units) || !any(selected_edges)) next

    counts <- table(edge_target[selected_edges])
    targets <- intersect(names(counts), rownames(degree))
    if (!length(targets)) next
    degree[targets, selected_units] <- matrix(
      as.numeric(counts[targets]),
      nrow = length(targets), ncol = length(selected_units)
    )
  }
  degree
}

.rc_average_sign_only_projection <- function(projection, degree) {
  projection <- as.matrix(projection)
  degree <- as.matrix(degree)
  if (!identical(dimnames(projection), dimnames(degree))) {
    stop(
      "Sign-only projection and target active-edge degree must align exactly.",
      call. = FALSE
    )
  }
  valid <- is.finite(projection) & is.finite(degree) & degree > 0
  projection[valid] <- projection[valid] / degree[valid]
  projection
}
