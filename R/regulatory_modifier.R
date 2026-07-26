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

.rc_regulatory_edge_projection_weight <- function(edges) {
  effective <- if ("effective_estimate" %in% colnames(edges)) {
    suppressWarnings(as.numeric(edges$effective_estimate))
  } else {
    suppressWarnings(as.numeric(edges$estimate))
  }
  stability <- if ("stability_weight" %in% colnames(edges)) {
    suppressWarnings(as.numeric(edges$stability_weight))
  } else {
    rep(1, nrow(edges))
  }
  tf_reference <- if ("tf_reference" %in% colnames(edges)) {
    suppressWarnings(as.numeric(edges$tf_reference))
  } else {
    rep(1, nrow(edges))
  }
  interaction_scale <- if ("interaction_scale" %in% colnames(edges)) {
    suppressWarnings(as.numeric(edges$interaction_scale))
  } else {
    rep(1, nrow(edges))
  }
  valid <- is.finite(effective) & is.finite(stability) & stability >= 0 &
    is.finite(tf_reference) & tf_reference >= 0 &
    is.finite(interaction_scale) & interaction_scale > 0
  answer <- rep(0, nrow(edges))
  answer[valid] <- effective[valid] * stability[valid] *
    tf_reference[valid] / interaction_scale[valid]
  answer
}

.rc_target_regulatory_reliability <- function(edges) {
  value <- if ("cv_rsq" %in% colnames(edges)) {
    suppressWarnings(as.numeric(edges$cv_rsq))
  } else if ("rsq" %in% colnames(edges)) {
    suppressWarnings(as.numeric(edges$rsq))
  } else {
    numeric()
  }
  value <- value[is.finite(value)]
  if (!length(value)) return(0)
  sqrt(min(max(stats::median(value), 0), 1))
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
      "Regulatory edge table is missing columns: ",
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
    "target reliability is sqrt(clipped median cross-validated R-squared)",
    "shared across conditions within one cell type"
  )
  attr(modifier, "projection_formula") <- paste(
    "edge weight = effective coefficient * stability * shared TF reference /",
    "interaction scale; TF edges sharing one ATAC peak are signed-summed;",
    "the denominator is max_condition sum_peak abs(weight)"
  )
  if (!nrow(significant_edges) || !length(genes)) return(modifier)

  edges <- significant_edges
  edges$target <- toupper(trimws(as.character(edges$target)))
  edges$tf <- toupper(trimws(as.character(edges$tf)))
  edges$region <- trimws(as.character(edges$region))
  edges$.projection_weight <- .rc_regulatory_edge_projection_weight(edges)
  edges <- edges[
    !is.na(edges$target) & nzchar(edges$target) &
      !is.na(edges$tf) & nzchar(edges$tf) &
      !is.na(edges$region) & nzchar(edges$region) &
      is.finite(edges$.projection_weight) & edges$.projection_weight != 0,
    , drop = FALSE
  ]
  if (!nrow(edges)) return(modifier)

  atac <- .rc_pando_assay_data(object, atac_assay)
  if ("atac_feature_id" %in% colnames(edges)) {
    supplied <- trimws(as.character(edges$atac_feature_id))
    edges$.peak_id <- ifelse(supplied %in% rownames(atac), supplied, NA_character_)
  } else {
    edges$.peak_id <- NA_character_
  }
  unresolved <- is.na(edges$.peak_id) | !nzchar(edges$.peak_id)
  if (any(unresolved)) {
    peak_keys <- toupper(.rc_pando_region_key(rownames(atac)))
    peak_keep <- !is.na(peak_keys) & nzchar(peak_keys) & !duplicated(peak_keys)
    peak_lookup <- stats::setNames(rownames(atac)[peak_keep], peak_keys[peak_keep])
    edges$.peak_id[unresolved] <- unname(
      peak_lookup[toupper(.rc_pando_region_key(edges$region[unresolved]))]
    )
  }
  edges <- edges[
    !is.na(edges$.peak_id) & nzchar(edges$.peak_id), , drop = FALSE
  ]
  if (!nrow(edges)) return(modifier)

  celltypes <- unique(as.character(edges[[celltype_col]]))
  for (cell_type in celltypes) {
    celltype_edges <- edges[
      as.character(edges[[celltype_col]]) == cell_type, , drop = FALSE
    ]
    reference_units <- units[
      as.character(unit_meta[[celltype_col]]) == cell_type
    ]
    if (!nrow(celltype_edges) || !length(reference_units)) next

    for (target in unique(celltype_edges$target)) {
      target_edges <- celltype_edges[
        celltype_edges$target == target, , drop = FALSE
      ]
      gene_id <- tolower(target)
      if (!gene_id %in% rownames(modifier) || !nrow(target_edges)) next
      reliability <- .rc_target_regulatory_reliability(target_edges)
      if (!is.finite(reliability) || reliability <= 0) next

      condition_peak_key <- paste(
        as.character(target_edges[[condition_col]]),
        target_edges$.peak_id,
        sep = "\001"
      )
      peak_rows <- split(seq_len(nrow(target_edges)), condition_peak_key)
      peak_weights <- do.call(rbind, lapply(peak_rows, function(index) {
        one <- target_edges[index, , drop = FALSE]
        data.frame(
          condition = as.character(one[[condition_col]][[1L]]),
          peak_id = as.character(one$.peak_id[[1L]]),
          weight = sum(one$.projection_weight),
          stringsAsFactors = FALSE
        )
      }))
      peak_weights <- peak_weights[
        is.finite(peak_weights$weight) & peak_weights$weight != 0,
        , drop = FALSE
      ]
      if (!nrow(peak_weights)) next
      denominator_by_condition <- vapply(
        split(abs(peak_weights$weight), peak_weights$condition),
        sum,
        numeric(1)
      )
      shared_denominator <- max(denominator_by_condition, na.rm = TRUE)
      if (!is.finite(shared_denominator) || shared_denominator <= 0) next

      peaks <- unique(peak_weights$peak_id)
      edge_activity <- rc_gene_score(
        as.matrix(atac[peaks, units, drop = FALSE]),
        mode = "absolute",
        half_saturation = getOption("RegCompassR.atac_half_saturation", 1)
      )
      rownames(edge_activity) <- peaks
      edge_deviation <- .rc_edge_activity_deviation(
        edge_activity[, reference_units, drop = FALSE],
        min_scale = min_scale
      )
      full_deviation <- matrix(
        NA_real_, nrow = length(peaks), ncol = length(units),
        dimnames = list(peaks, units)
      )
      full_deviation[, reference_units] <- edge_deviation

      for (condition in unique(peak_weights$condition)) {
        group_units <- units[
          as.character(unit_meta[[celltype_col]]) == cell_type &
            as.character(unit_meta[[condition_col]]) == condition
        ]
        if (!length(group_units)) next
        one <- peak_weights[
          peak_weights$condition == condition, , drop = FALSE
        ]
        weights <- stats::setNames(
          one$weight / shared_denominator, one$peak_id
        )
        value <- reliability * as.numeric(crossprod(
          weights,
          full_deviation[names(weights), group_units, drop = FALSE]
        ))
        modifier[gene_id, group_units] <- pmax(pmin(value, 1), -1)
      }
    }
  }
  modifier
}
