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

.rc_case_insensitive_lookup <- function(ids) {
  keys <- toupper(trimws(as.character(ids)))
  keep <- !is.na(keys) & nzchar(keys) & !duplicated(keys)
  stats::setNames(as.character(ids)[keep], keys[keep])
}

.rc_pando_gene_confidence <- function(significant_edges, object, atac_assay,
                                      target_genes = NULL,
                                      rna_assay = "RNA") {
  units <- colnames(object)
  genes <- unique(toupper(as.character(target_genes)))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (!length(genes) && is.data.frame(significant_edges) &&
      "target" %in% colnames(significant_edges)) {
    genes <- unique(toupper(as.character(significant_edges$target)))
    genes <- genes[!is.na(genes) & nzchar(genes)]
  }
  confidence <- matrix(
    NA_real_, nrow = length(genes), ncol = length(units),
    dimnames = list(tolower(genes), units)
  )
  diagnostics <- data.frame(
    gene = genes,
    n_pando_edges = integer(length(genes)),
    n_positive_edges = integer(length(genes)),
    n_negative_edges = integer(length(genes)),
    n_unique_regions = integer(length(genes)),
    n_matched_regions = integer(length(genes)),
    n_unique_tfs = integer(length(genes)),
    n_matched_tfs = integer(length(genes)),
    matched_region_fraction = NA_real_,
    matched_tf_fraction = NA_real_,
    pando_supported = FALSE,
    confidence_source = "pando_signed_tf_peak_gene_regulatory_support",
    stringsAsFactors = FALSE
  )
  if (!is.data.frame(significant_edges) || !nrow(significant_edges)) {
    return(list(gene_confidence = confidence, diagnostics = diagnostics))
  }
  required <- c("target", "region", "tf", "estimate")
  missing <- setdiff(required, colnames(significant_edges))
  if (length(missing)) {
    stop("Pando coefficient table lacks columns required for signed support: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  edges <- significant_edges
  edges$target <- toupper(trimws(as.character(edges$target)))
  edges$tf <- toupper(trimws(as.character(edges$tf)))
  edges$region <- trimws(as.character(edges$region))
  edges$estimate <- suppressWarnings(as.numeric(edges$estimate))
  edges <- edges[
    !is.na(edges$target) & nzchar(edges$target) &
      !is.na(edges$tf) & nzchar(edges$tf) &
      !is.na(edges$region) & nzchar(edges$region) &
      is.finite(edges$estimate) & edges$estimate != 0 &
      edges$target %in% genes, , drop = FALSE
  ]
  if (!nrow(edges)) {
    return(list(gene_confidence = confidence, diagnostics = diagnostics))
  }
  edges$.sign <- sign(edges$estimate)
  edges$.weight <- abs(edges$estimate)
  if ("rsq" %in% colnames(edges)) {
    rsq <- suppressWarnings(as.numeric(edges$rsq))
    quality <- sqrt(pmax(rsq, 0))
    quality[!is.finite(quality)] <- 0
    edges$.weight <- edges$.weight * quality
  }
  edges <- edges[is.finite(edges$.weight) & edges$.weight > 0, , drop = FALSE]
  if (!nrow(edges)) {
    return(list(gene_confidence = confidence, diagnostics = diagnostics))
  }
  atac <- .rc_pando_assay_data(object, atac_assay)
  rna <- .rc_pando_assay_data(object, rna_assay)
  peak_keys <- toupper(.rc_pando_region_key(rownames(atac)))
  keep_peak <- !is.na(peak_keys) & nzchar(peak_keys) & !duplicated(peak_keys)
  peak_lookup <- stats::setNames(rownames(atac)[keep_peak], peak_keys[keep_peak])
  edges$.peak_id <- unname(
    peak_lookup[toupper(.rc_pando_region_key(edges$region))]
  )
  tf_lookup <- .rc_case_insensitive_lookup(rownames(rna))
  edges$.tf_id <- unname(tf_lookup[edges$tf])

  for (gene in genes) {
    selected <- edges[edges$target == gene, , drop = FALSE]
    row <- diagnostics$gene == gene
    diagnostics$n_pando_edges[row] <- nrow(selected)
    diagnostics$n_positive_edges[row] <- sum(selected$.sign > 0)
    diagnostics$n_negative_edges[row] <- sum(selected$.sign < 0)
    diagnostics$n_unique_regions[row] <- length(unique(selected$region))
    diagnostics$n_unique_tfs[row] <- length(unique(selected$tf))
    matched <- selected[
      !is.na(selected$.peak_id) & nzchar(selected$.peak_id) &
        !is.na(selected$.tf_id) & nzchar(selected$.tf_id), , drop = FALSE
    ]
    diagnostics$n_matched_regions[row] <- length(unique(matched$.peak_id))
    diagnostics$n_matched_tfs[row] <- length(unique(matched$.tf_id))
    diagnostics$matched_region_fraction[row] <- if (
      diagnostics$n_unique_regions[row] > 0L
    ) diagnostics$n_matched_regions[row] / diagnostics$n_unique_regions[row] else
      NA_real_
    diagnostics$matched_tf_fraction[row] <- if (
      diagnostics$n_unique_tfs[row] > 0L
    ) diagnostics$n_matched_tfs[row] / diagnostics$n_unique_tfs[row] else
      NA_real_
    diagnostics$pando_supported[row] <- nrow(matched) > 0L
    if (!nrow(matched)) next
    peak_score <- rc_gene_score(
      as.matrix(atac[matched$.peak_id, units, drop = FALSE]),
      mode = "absolute",
      half_saturation = getOption("RegCompassR.atac_half_saturation", 1)
    )
    tf_score <- rc_gene_score(
      as.matrix(rna[matched$.tf_id, units, drop = FALSE]),
      mode = "absolute",
      half_saturation = getOption("RegCompassR.tf_half_saturation", 1)
    )
    edge_activity <- peak_score * tf_score
    weights <- matched$.weight / sum(matched$.weight)
    signed_activity <- as.numeric(crossprod(
      weights * matched$.sign, edge_activity
    ))
    confidence[tolower(gene), ] <- .rc_clamp01(
      0.5 + 0.5 * signed_activity
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
                                           atac_assay = "ATAC",
                                           rna_assay = "RNA") {
  parsed <- rc_parse_gpr_table(gem$gpr_table)
  target_genes <- meta_modules$target_metabolic_genes %||%
    rc_metabolic_gpr_genes(gem$gpr_table)
  gene <- .rc_pando_gene_confidence(
    meta_modules$tf_peak_gene_significant,
    pando_object,
    atac_assay = atac_assay,
    rna_assay = rna_assay,
    target_genes = target_genes
  )
  supported <- rowSums(is.finite(gene$gene_confidence)) > 0L
  gene_for_reaction <- gene$gene_confidence[supported, , drop = FALSE]
  reaction <- rc_reaction_confidence(
    parsed,
    gene_confidence = gene_for_reaction,
    and_method = "min",
    or_method = "max",
    unit_ids = colnames(pando_object)
  )
  list(
    gene_confidence = gene$gene_confidence,
    gene_confidence_diagnostics = gene$diagnostics,
    reaction_confidence = reaction,
    reaction_confidence_matrix = .rc_pando_reaction_confidence_matrix(
      reaction, names(parsed), colnames(pando_object)
    ),
    confidence_source = "pando_signed_tf_peak_gene_regulatory_support"
  )
}
