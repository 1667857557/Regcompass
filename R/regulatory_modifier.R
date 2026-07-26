.rc_edge_activity_deviation <- function(edge_activity, min_scale = 0.05) {
  edge_activity <- as.matrix(edge_activity)
  centers <- matrixStats::rowMedians(edge_activity, na.rm = TRUE)
  mad_scale <- matrixStats::rowMads(
    edge_activity, constant = 1.4826, na.rm = TRUE
  )
  iqr_scale <- matrixStats::rowIQRs(edge_activity, na.rm = TRUE) / 1.349
  scale <- pmax(mad_scale, iqr_scale, min_scale, na.rm = TRUE)
  standardized <- sweep(edge_activity, 1L, centers, "-")
  standardized <- sweep(standardized, 1L, scale, "/")
  tanh(standardized)
}

.rc_condition_gene_regulatory_modifier <- function(
    significant_edges, object, unit_meta,
    condition_col = "condition", celltype_col = "cell_type",
    atac_assay = "ATAC", target_genes = NULL, min_scale = 0.05) {
  if (!is.data.frame(significant_edges)) {
    stop("`significant_edges` must be a data.frame.", call. = FALSE)
  }
  required_edges <- c(
    "target", "region", "tf", "estimate", condition_col, celltype_col
  )
  missing_edges <- setdiff(required_edges, colnames(significant_edges))
  if (length(missing_edges)) {
    stop(
      "GRN edge table is missing columns: ",
      paste(missing_edges, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.data.frame(unit_meta) ||
      !all(c("pool_id", condition_col, celltype_col) %in% colnames(unit_meta))) {
    stop(
      "`unit_meta` is incomplete for condition-pooled regulatory scoring.",
      call. = FALSE
    )
  }
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
    "joint cross-validated target R-squared is preferred; legacy finite",
    "Pando R-squared is used only when cross-validated R-squared is absent"
  )
  attr(modifier, "projection_policy") <- paste(
    "TF-peak-target stable coefficients are converted to ATAC-only projection",
    "weights; TFs sharing one peak are signed-summed; each target uses one",
    "cell-type-wide denominator shared across conditions"
  )
  if (!nrow(significant_edges) || !length(genes)) return(modifier)

  edges <- significant_edges
  edges$target <- toupper(trimws(as.character(edges$target)))
  edges$tf <- toupper(trimws(as.character(edges$tf)))
  edges$region <- trimws(as.character(edges$region))
  edges$estimate <- suppressWarnings(as.numeric(edges$estimate))
  projection_weight <- if ("atac_projection_weight" %in% colnames(edges)) {
    suppressWarnings(as.numeric(edges$atac_projection_weight))
  } else {
    edges$estimate
  }
  edges$.projection_weight <- projection_weight
  edges <- edges[
    !is.na(edges$target) & nzchar(edges$target) &
      !is.na(edges$tf) & nzchar(edges$tf) &
      !is.na(edges$region) & nzchar(edges$region) &
      is.finite(edges$estimate) & is.finite(edges$.projection_weight) &
      edges$.projection_weight != 0,
    , drop = FALSE
  ]
  if (!nrow(edges)) return(modifier)

  atac <- .rc_pando_assay_data(object, atac_assay)
  if ("atac_feature_id" %in% colnames(edges)) {
    feature <- trimws(as.character(edges$atac_feature_id))
    feature[!feature %in% rownames(atac)] <- NA_character_
    edges$.peak_id <- feature
  } else {
    peak_keys <- toupper(.rc_pando_region_key(rownames(atac)))
    peak_keep <- !is.na(peak_keys) & nzchar(peak_keys) & !duplicated(peak_keys)
    peak_lookup <- stats::setNames(rownames(atac)[peak_keep], peak_keys[peak_keep])
    edges$.peak_id <- unname(
      peak_lookup[toupper(.rc_pando_region_key(edges$region))]
    )
  }
  edges <- edges[
    !is.na(edges$.peak_id) & nzchar(edges$.peak_id), , drop = FALSE
  ]
  if (!nrow(edges)) return(modifier)

  reliability_column <- if ("cv_rsq" %in% colnames(edges)) {
    "cv_rsq"
  } else if ("rsq" %in% colnames(edges)) {
    "rsq"
  } else {
    NULL
  }
  celltypes <- unique(as.character(edges[[celltype_col]]))
  for (celltype_value in celltypes) {
    celltype_edges <- edges[
      as.character(edges[[celltype_col]]) == celltype_value,
      , drop = FALSE
    ]
    reference_units <- units[
      as.character(unit_meta[[celltype_col]]) == celltype_value
    ]
    if (!nrow(celltype_edges) || !length(reference_units)) next

    for (target in unique(celltype_edges$target)) {
      target_edges <- celltype_edges[
        celltype_edges$target == target, , drop = FALSE
      ]
      gene_id <- tolower(target)
      if (!gene_id %in% rownames(modifier) || !nrow(target_edges)) next
      condition_edges <- split(
        target_edges,
        as.character(target_edges[[condition_col]])
      )
      collapsed <- lapply(condition_edges, function(one) {
        value <- tapply(
          one$.projection_weight,
          one$.peak_id,
          sum,
          na.rm = TRUE
        )
        value <- value[is.finite(value) & value != 0]
        value
      })
      denominator <- max(vapply(collapsed, function(value) {
        sum(abs(value))
      }, numeric(1)), na.rm = TRUE)
      if (!is.finite(denominator) || denominator <= 0) next
      peak_ids <- unique(unlist(lapply(collapsed, names), use.names = FALSE))
      peak_ids <- intersect(peak_ids, rownames(atac))
      if (!length(peak_ids)) next

      edge_activity <- rc_gene_score(
        as.matrix(atac[peak_ids, units, drop = FALSE]),
        mode = "absolute",
        half_saturation = getOption("RegCompassR.atac_half_saturation", 1)
      )
      edge_deviation_reference <- .rc_edge_activity_deviation(
        edge_activity[, reference_units, drop = FALSE],
        min_scale = min_scale
      )
      reliability_values <- if (!is.null(reliability_column)) {
        suppressWarnings(as.numeric(target_edges[[reliability_column]]))
      } else {
        numeric()
      }
      reliability_values <- reliability_values[is.finite(reliability_values)]
      reliability <- if (length(reliability_values)) {
        sqrt(min(max(stats::median(reliability_values), 0), 1))
      } else {
        0
      }

      for (condition_value in names(collapsed)) {
        group_units <- units[
          as.character(unit_meta[[celltype_col]]) == celltype_value &
            as.character(unit_meta[[condition_col]]) == condition_value
        ]
        if (!length(group_units)) next
        weight <- stats::setNames(rep(0, length(peak_ids)), peak_ids)
        value <- collapsed[[condition_value]]
        common <- intersect(names(value), peak_ids)
        weight[common] <- value[common] / denominator
        edge_deviation <- edge_deviation_reference[
          peak_ids, group_units, drop = FALSE
        ]
        score <- reliability * as.numeric(crossprod(weight, edge_deviation))
        modifier[gene_id, group_units] <- pmax(pmin(score, 1), -1)
      }
    }
  }
  modifier
}
