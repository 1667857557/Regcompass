# Dispatch regulatory projection by GRN evidence schema.

.rc_condition_gene_regulatory_modifier_shared <-
  .rc_condition_gene_regulatory_modifier

.rc_condition_gene_regulatory_modifier_legacy <- function(
    significant_edges, object, unit_meta,
    condition_col = "condition", celltype_col = "cell_type",
    atac_assay = "ATAC", target_genes = NULL, min_scale = 0.05) {
  units <- colnames(object)
  unit_meta <- unit_meta[
    match(units, as.character(unit_meta$pool_id)), , drop = FALSE
  ]
  if (anyNA(unit_meta$pool_id)) {
    stop("Metacell metadata do not align to the scoring object.", call. = FALSE)
  }
  genes <- unique(tolower(trimws(as.character(target_genes))))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (!length(genes)) {
    genes <- unique(tolower(trimws(as.character(significant_edges$target))))
  }
  modifier <- matrix(
    0,
    nrow = length(genes),
    ncol = length(units),
    dimnames = list(genes, units)
  )
  attr(modifier, "reliability_policy") <- paste(
    "legacy condition-Pando projection: finite Pando R-squared values only;",
    "targets without finite R-squared receive reliability zero"
  )
  attr(modifier, "projection_formula") <- paste(
    "legacy per-condition target weights equal abs(Pando estimate) divided by",
    "their condition-target sum, multiplied by coefficient sign"
  )
  if (!nrow(significant_edges) || !length(genes)) return(modifier)

  edges <- significant_edges
  rsq <- if ("rsq" %in% colnames(edges)) {
    suppressWarnings(as.numeric(edges$rsq))
  } else {
    rep(NA_real_, nrow(edges))
  }
  edges <- edges[is.finite(rsq), , drop = FALSE]
  if (!nrow(edges)) return(modifier)
  edges$target <- toupper(trimws(as.character(edges$target)))
  edges$tf <- toupper(trimws(as.character(edges$tf)))
  edges$region <- trimws(as.character(edges$region))
  edges$estimate <- suppressWarnings(as.numeric(edges$estimate))
  edges <- edges[
    !is.na(edges$target) & nzchar(edges$target) &
      !is.na(edges$tf) & nzchar(edges$tf) &
      !is.na(edges$region) & nzchar(edges$region) &
      is.finite(edges$estimate) & edges$estimate != 0,
    , drop = FALSE
  ]
  if (!nrow(edges)) return(modifier)

  atac <- .rc_pando_assay_data(object, atac_assay)
  peak_keys <- toupper(.rc_pando_region_key(rownames(atac)))
  peak_keep <- !is.na(peak_keys) & nzchar(peak_keys) & !duplicated(peak_keys)
  peak_lookup <- stats::setNames(rownames(atac)[peak_keep], peak_keys[peak_keep])
  edges$.peak_id <- unname(
    peak_lookup[toupper(.rc_pando_region_key(edges$region))]
  )
  edges <- edges[
    !is.na(edges$.peak_id) & nzchar(edges$.peak_id), , drop = FALSE
  ]
  if (!nrow(edges)) return(modifier)

  group_key_edges <- paste(
    as.character(edges[[condition_col]]),
    as.character(edges[[celltype_col]]),
    sep = "\001"
  )
  group_key_units <- paste(
    as.character(unit_meta[[condition_col]]),
    as.character(unit_meta[[celltype_col]]),
    sep = "\001"
  )
  for (group_key in unique(group_key_edges)) {
    group_edges <- edges[group_key_edges == group_key, , drop = FALSE]
    group_units <- units[group_key_units == group_key]
    if (!nrow(group_edges) || !length(group_units)) next
    for (target in unique(group_edges$target)) {
      selected <- group_edges[group_edges$target == target, , drop = FALSE]
      gene_id <- tolower(target)
      if (!gene_id %in% rownames(modifier) || !nrow(selected)) next
      edge_activity <- rc_gene_score(
        as.matrix(atac[selected$.peak_id, units, drop = FALSE]),
        mode = "absolute",
        half_saturation = getOption("RegCompassR.atac_half_saturation", 1)
      )
      celltype_value <- as.character(group_edges[[celltype_col]][[1L]])
      reference_units <- units[
        as.character(unit_meta[[celltype_col]]) == celltype_value
      ]
      if (!length(reference_units)) next
      edge_deviation_reference <- .rc_edge_activity_deviation(
        edge_activity[, reference_units, drop = FALSE],
        min_scale = min_scale
      )
      edge_deviation <- edge_deviation_reference[, group_units, drop = FALSE]
      weight <- abs(selected$estimate)
      weight[!is.finite(weight)] <- 0
      if (!any(weight > 0)) next
      weight <- weight / sum(weight)
      model_rsq <- suppressWarnings(as.numeric(selected$rsq))
      model_rsq <- model_rsq[is.finite(model_rsq)]
      reliability <- if (length(model_rsq)) {
        sqrt(min(max(stats::median(model_rsq), 0), 1))
      } else {
        0
      }
      signed_weight <- weight * sign(selected$estimate)
      value <- reliability * as.numeric(
        crossprod(signed_weight, edge_deviation)
      )
      modifier[gene_id, group_units] <- pmax(pmin(value, 1), -1)
    }
  }
  modifier
}

.rc_condition_gene_regulatory_modifier <- function(
    significant_edges, object, unit_meta,
    condition_col = "condition", celltype_col = "cell_type",
    atac_assay = "ATAC", target_genes = NULL, min_scale = 0.05) {
  multitask_columns <- c(
    "effective_estimate", "interaction_scale", "tf_reference",
    "stability_weight"
  )
  shared_mode <- is.data.frame(significant_edges) &&
    all(multitask_columns %in% colnames(significant_edges))
  function_to_call <- if (shared_mode) {
    .rc_condition_gene_regulatory_modifier_shared
  } else {
    .rc_condition_gene_regulatory_modifier_legacy
  }
  function_to_call(
    significant_edges = significant_edges,
    object = object,
    unit_meta = unit_meta,
    condition_col = condition_col,
    celltype_col = celltype_col,
    atac_assay = atac_assay,
    target_genes = target_genes,
    min_scale = min_scale
  )
}
