.rc_pando_assay_data <- function(object, assay) {
  value <- tryCatch(
    SeuratObject::GetAssayData(object, assay = assay, slot = "data"),
    error = function(e) NULL
  )
  if (is.null(value)) {
    value <- tryCatch(
      SeuratObject::GetAssayData(object, assay = assay, layer = "data"),
      error = function(e) NULL
    )
  }
  if (is.null(value) || nrow(value) == 0L || ncol(value) == 0L) {
    stop(
      "Pando-derived confidence requires normalized ATAC data in the Pando metacell object.",
      call. = FALSE
    )
  }
  value
}

.rc_pando_region_key <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^([^:]+):(\\d+)-(\\d+)$", "\\1-\\2-\\3", x)
  x <- sub("^([^:]+):(\\d+):(\\d+)$", "\\1-\\2-\\3", x)
  x
}

.rc_pando_gene_confidence <- function(significant_edges, object, atac_assay,
                                      target_genes = NULL) {
  units <- colnames(object)
  genes <- unique(toupper(as.character(target_genes)))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (!length(genes) && is.data.frame(significant_edges) &&
      "target" %in% colnames(significant_edges)) {
    genes <- unique(toupper(as.character(significant_edges$target)))
    genes <- genes[!is.na(genes) & nzchar(genes)]
  }
  confidence <- matrix(
    NA_real_,
    nrow = length(genes),
    ncol = length(units),
    dimnames = list(tolower(genes), units)
  )
  diagnostics <- data.frame(
    gene = genes,
    n_pando_edges = integer(length(genes)),
    n_unique_regions = integer(length(genes)),
    n_matched_regions = integer(length(genes)),
    matched_region_fraction = NA_real_,
    pando_supported = FALSE,
    confidence_source = "pando_internal_peak_gene_accessibility",
    stringsAsFactors = FALSE
  )
  if (!is.data.frame(significant_edges) || !nrow(significant_edges)) {
    return(list(gene_confidence = confidence, diagnostics = diagnostics))
  }
  required <- c("target", "region")
  missing <- setdiff(required, colnames(significant_edges))
  if (length(missing)) {
    stop(
      "Pando coefficient table lacks columns required for confidence: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  edges <- significant_edges
  edges$target <- toupper(trimws(as.character(edges$target)))
  edges$region <- trimws(as.character(edges$region))
  edges <- edges[
    !is.na(edges$target) & nzchar(edges$target) &
      !is.na(edges$region) & nzchar(edges$region) &
      edges$target %in% genes,
    ,
    drop = FALSE
  ]
  if (!nrow(edges)) {
    return(list(gene_confidence = confidence, diagnostics = diagnostics))
  }
  edge_weight <- rep(1, nrow(edges))
  if ("estimate" %in% colnames(edges)) {
    estimate <- abs(suppressWarnings(as.numeric(edges$estimate)))
    use <- is.finite(estimate) & estimate > 0
    edge_weight[use] <- estimate[use]
  }
  if ("rsq" %in% colnames(edges)) {
    rsq <- suppressWarnings(as.numeric(edges$rsq))
    quality <- sqrt(pmax(rsq, 0))
    quality[!is.finite(quality) | quality <= 0] <- 1
    edge_weight <- edge_weight * quality
  }
  edge_weight[!is.finite(edge_weight) | edge_weight <= 0] <- 1
  edges$.weight <- edge_weight

  atac <- .rc_pando_assay_data(object, atac_assay)
  peak_ids <- rownames(atac)
  peak_keys <- .rc_pando_region_key(peak_ids)
  keep_peak <- !duplicated(peak_keys)
  peak_lookup <- stats::setNames(peak_ids[keep_peak], peak_keys[keep_peak])
  edges$.peak_id <- unname(peak_lookup[.rc_pando_region_key(edges$region)])

  for (gene in genes) {
    selected <- edges[edges$target == gene, , drop = FALSE]
    row <- diagnostics$gene == gene
    diagnostics$n_pando_edges[row] <- nrow(selected)
    diagnostics$n_unique_regions[row] <- length(unique(selected$region))
    matched <- selected[
      !is.na(selected$.peak_id) & nzchar(selected$.peak_id),
      ,
      drop = FALSE
    ]
    diagnostics$n_matched_regions[row] <- length(unique(matched$.peak_id))
    denominator <- diagnostics$n_unique_regions[row]
    diagnostics$matched_region_fraction[row] <- if (denominator > 0L) {
      diagnostics$n_matched_regions[row] / denominator
    } else {
      NA_real_
    }
    diagnostics$pando_supported[row] <- nrow(matched) > 0L
    if (!nrow(matched)) next
    aggregated <- stats::aggregate(
      matched$.weight,
      by = list(peak_id = matched$.peak_id),
      FUN = sum,
      na.rm = TRUE
    )
    colnames(aggregated)[[2L]] <- "weight"
    accessibility <- as.matrix(
      atac[aggregated$peak_id, units, drop = FALSE]
    )
    accessibility_score <- rc_gene_score(accessibility)
    accessibility_score[!is.finite(accessibility) | accessibility <= 0] <- 0
    weights <- aggregated$weight
    weights[!is.finite(weights) | weights <= 0] <- 1
    weights <- weights / sum(weights)
    confidence[tolower(gene), ] <- as.numeric(
      crossprod(weights, accessibility_score)
    )
  }
  list(gene_confidence = confidence, diagnostics = diagnostics)
}

.rc_pando_reaction_confidence_matrix <- function(confidence, reaction_ids,
                                                 unit_ids) {
  out <- matrix(
    NA_real_,
    nrow = length(reaction_ids),
    ncol = length(unit_ids),
    dimnames = list(reaction_ids, unit_ids)
  )
  if (!is.data.frame(confidence) || !nrow(confidence)) return(out)
  required <- c("reaction_id", "pool_id", "reaction_confidence")
  if (!all(required %in% colnames(confidence))) {
    stop("Pando reaction confidence table is incomplete.", call. = FALSE)
  }
  row_index <- match(as.character(confidence$reaction_id), reaction_ids)
  col_index <- match(as.character(confidence$pool_id), unit_ids)
  keep <- !is.na(row_index) & !is.na(col_index)
  out[cbind(row_index[keep], col_index[keep])] <- as.numeric(
    confidence$reaction_confidence[keep]
  )
  out
}

.rc_pando_reaction_confidence <- function(meta_modules, pando_object, gem,
                                          atac_assay = "ATAC") {
  parsed <- rc_parse_gpr_table(gem$gpr_table)
  target_genes <- meta_modules$target_metabolic_genes %||%
    rc_metabolic_gpr_genes(gem$gpr_table)
  gene <- .rc_pando_gene_confidence(
    meta_modules$tf_peak_gene_significant,
    pando_object,
    atac_assay = atac_assay,
    target_genes = target_genes
  )
  supported <- rowSums(is.finite(gene$gene_confidence)) > 0L
  gene_for_reaction <- gene$gene_confidence[
    supported,
    ,
    drop = FALSE
  ]
  reaction <- rc_reaction_confidence(
    parsed,
    gene_confidence = gene_for_reaction,
    unit_ids = colnames(pando_object)
  )
  list(
    gene_confidence = gene$gene_confidence,
    gene_confidence_diagnostics = gene$diagnostics,
    reaction_confidence = reaction,
    reaction_confidence_matrix = .rc_pando_reaction_confidence_matrix(
      reaction,
      names(parsed),
      colnames(pando_object)
    ),
    confidence_source = "pando_internal_peak_gene_accessibility"
  )
}
