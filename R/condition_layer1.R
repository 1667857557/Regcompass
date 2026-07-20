.rc_integrate_regulatory_support_v170 <- function(
    rna_support, regulatory_modifier, alpha = 1) {
  rna_support <- as.matrix(rna_support)
  regulatory_modifier <- as.matrix(regulatory_modifier)
  if (!identical(dim(rna_support), dim(regulatory_modifier)) ||
      !identical(dimnames(rna_support), dimnames(regulatory_modifier))) {
    stop("RNA support and regulatory modifier matrices must align exactly.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha < 0) {
    stop("`alpha` must be one finite non-negative number.", call. = FALSE)
  }
  C <- pmin(pmax(rna_support, 0), 1)
  R <- pmin(pmax(regulatory_modifier, -1), 1)
  multiplier <- 2^(alpha * R)
  numerator <- C * multiplier
  denominator <- 1 - C + numerator
  out <- numerator / denominator
  out[C <= 0] <- 0
  out[C >= 1] <- 1
  out[!is.finite(out)] <- NA_real_
  dimnames(out) <- dimnames(C)
  attr(out, "integration_formula") <- paste(
    "C_multiome = C_RNA * 2^(alpha * R_ATAC) /",
    "(1 - C_RNA + C_RNA * 2^(alpha * R_ATAC))"
  )
  attr(out, "score_semantics") <- paste(
    "zero-preserving bounded gene support with a signed accessibility-derived",
    "regulatory modifier on the support log-odds scale"
  )
  out
}

.rc_build_condition_pooled_layer1 <- function(
    metacell_object, meta_modules, gem, metacell_meta,
    sample_col = "sample_id", condition_col = "condition",
    celltype_col = "cell_type", rna_assay = "RNA", atac_assay = "ATAC",
    regulatory_alpha = 1, gpr_tau = 0.20,
    gene_half_saturation = getOption("RegCompassR.cpm_half_saturation", 1)) {
  parsed <- rc_parse_gpr_table(gem$gpr_table)
  gpr_genes <- unique(tolower(unlist(parsed, use.names = FALSE)))
  counts <- .rc_get_assay_counts(metacell_object, rna_assay)
  full_library_size <- Matrix::colSums(counts)
  keep <- tolower(rownames(counts)) %in% gpr_genes
  rna_counts <- counts[keep, , drop = FALSE]
  rna_logcpm <- .rc_metacell_logcpm(
    rna_counts,
    library_size = full_library_size[colnames(rna_counts)]
  )
  rownames(rna_logcpm) <- tolower(rownames(rna_logcpm))
  if (anyDuplicated(rownames(rna_logcpm))) {
    stop("Duplicated GPR gene identifiers after case normalization.", call. = FALSE)
  }

  unit_meta <- metacell_meta
  id_col <- if ("metacell_id" %in% colnames(unit_meta)) {
    "metacell_id"
  } else if ("pool_id" %in% colnames(unit_meta)) {
    "pool_id"
  } else {
    stop("Pooled metacell metadata lack metacell_id/pool_id.", call. = FALSE)
  }
  unit_meta$pool_id <- as.character(unit_meta[[id_col]])
  unit_meta$unit_id <- unit_meta$pool_id
  unit_meta[[sample_col]] <- paste0(
    as.character(unit_meta[[condition_col]]),
    "__pooled"
  )
  unit_meta <- unit_meta[
    match(colnames(rna_logcpm), unit_meta$pool_id),
    , drop = FALSE
  ]
  if (anyNA(unit_meta$pool_id)) {
    stop("Pooled metacell metadata do not align with RNA counts.", call. = FALSE)
  }

  gene_rna_support <- rc_gene_score(
    rna_logcpm,
    mode = "absolute",
    half_saturation = gene_half_saturation
  )
  modifier <- .rc_condition_gene_regulatory_modifier(
    significant_edges = meta_modules$tf_peak_gene_significant,
    object = metacell_object,
    unit_meta = unit_meta,
    condition_col = condition_col,
    celltype_col = celltype_col,
    atac_assay = atac_assay,
    target_genes = rownames(gene_rna_support)
  )
  modifier <- modifier[
    rownames(gene_rna_support),
    colnames(gene_rna_support),
    drop = FALSE
  ]
  gene_multiome_support <- .rc_integrate_regulatory_support_v170(
    gene_rna_support,
    modifier,
    alpha = regulatory_alpha
  )
  reaction_expression <- rc_reaction_capacity(
    parsed,
    gene_multiome_support,
    promiscuity_mode = "none",
    tau = gpr_tau,
    and_method = "boltzmann",
    or_method = "sum",
    BPPARAM = FALSE
  )

  list(
    schema_version = "regcompass_condition_pooled_layer1_v1.7.0",
    reaction_expression = reaction_expression,
    rna_metacell_logcpm = rna_logcpm,
    gene_support_rna = gene_rna_support,
    gene_regulatory_modifier = modifier,
    gene_support_multiome = gene_multiome_support,
    parsed_gpr = parsed,
    gpr_diagnostics = rc_gpr_diagnostics(parsed, rownames(rna_logcpm)),
    unit_meta = unit_meta,
    metacell_meta = unit_meta,
    layer1_unit = "condition_pooled_metacell",
    capacity_params = list(
      regulatory_alpha = regulatory_alpha,
      gene_half_saturation = gene_half_saturation,
      promiscuity_mode = "none",
      and_method = "boltzmann",
      tau = gpr_tau,
      or_method = "sum"
    ),
    evidence_formula = paste(
      "Pando coefficient-weighted ATAC accessibility modifier ->",
      "zero-preserving RNA support log-odds update ->",
      "Boltzmann AND and additive isozyme OR"
    ),
    evidence_provenance = list(
      direct_metacell_evidence = c("target_gene_RNA", "peak_ATAC"),
      learned_parameters = "condition_x_celltype Pando coefficients fit from RNA+ATAC",
      excluded_duplicate_evidence = "metacell TF RNA is not multiplied into the regulatory modifier",
      circularity_scope = paste(
        "Pando parameters remain estimated from the same pooled dataset; outputs",
        "are descriptive unless external fitting or cross-fitting is supplied"
      )
    )
  )
}

.rc_condition_penalty_comparison <- function(
    microcompass, condition_col = "condition", celltype_col = "cell_type",
    eps = 1e-8) {
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("`eps` must be one positive finite number.", call. = FALSE)
  }
  penalty <- as.matrix(microcompass$penalty)
  meta <- microcompass$unit_meta
  required <- c(condition_col, celltype_col)
  if (!is.numeric(penalty) || is.null(rownames(penalty)) ||
      is.null(colnames(penalty)) || anyDuplicated(rownames(penalty)) ||
      anyDuplicated(colnames(penalty))) {
    stop("microCOMPASS penalty must be a numeric matrix with unique dimnames.",
         call. = FALSE)
  }
  if (!is.data.frame(meta) || !all(required %in% colnames(meta))) {
    stop("microCOMPASS unit metadata lack condition/cell-type columns.", call. = FALSE)
  }
  unit_id <- if ("unit_id" %in% colnames(meta)) {
    as.character(meta$unit_id)
  } else if ("pool_id" %in% colnames(meta)) {
    as.character(meta$pool_id)
  } else {
    stop("microCOMPASS unit metadata lack unit_id/pool_id.", call. = FALSE)
  }
  if (anyNA(unit_id) || any(!nzchar(trimws(unit_id))) || anyDuplicated(unit_id)) {
    stop("microCOMPASS unit IDs must be unique and non-empty.", call. = FALSE)
  }
  if (!setequal(colnames(penalty), unit_id)) {
    stop("microCOMPASS penalties and unit metadata contain different units.",
         call. = FALSE)
  }
  meta$unit_id <- unit_id
  meta <- meta[match(colnames(penalty), meta$unit_id), , drop = FALSE]
  if (anyNA(meta[[condition_col]]) || anyNA(meta[[celltype_col]]) ||
      any(!nzchar(trimws(as.character(meta[[condition_col]])))) ||
      any(!nzchar(trimws(as.character(meta[[celltype_col]]))))) {
    stop("microCOMPASS condition/cell-type metadata are incomplete.", call. = FALSE)
  }

  row_meta <- rc_parse_microcompass_row_id(rownames(penalty))
  row_meta$row_id <- rownames(penalty)
  strata <- unique(meta[, c(condition_col, celltype_col), drop = FALSE])
  summary_rows <- lapply(seq_len(nrow(strata)), function(i) {
    condition <- as.character(strata[[condition_col]][[i]])
    celltype <- as.character(strata[[celltype_col]][[i]])
    keep <- as.character(meta[[condition_col]]) == condition &
      as.character(meta[[celltype_col]]) == celltype
    med <- matrixStats::rowMedians(penalty[, keep, drop = FALSE], na.rm = TRUE)
    med[is.nan(med)] <- NA_real_
    priority_rank <- rank(med, ties.method = "min", na.last = "keep")
    data.frame(
      row_id = row_meta$row_id,
      reaction_id = row_meta$reaction_id,
      target_direction = row_meta$target_direction,
      medium_scenario = row_meta$medium_scenario,
      condition = condition,
      cell_type = celltype,
      median_penalty = med,
      support_score = -log(med + eps),
      priority_rank = as.integer(priority_rank),
      n_metacells = sum(keep),
      descriptive_only = TRUE,
      biological_replicate_inference = FALSE,
      stringsAsFactors = FALSE
    )
  })
  ranking <- do.call(rbind, summary_rows)
  ranking <- ranking[order(
    ranking$cell_type,
    ranking$condition,
    ranking$priority_rank,
    ranking$reaction_id,
    ranking$target_direction,
    ranking$medium_scenario,
    na.last = TRUE
  ), , drop = FALSE]
  rownames(ranking) <- NULL

  contrast_rows <- list()
  contrast_index <- 1L
  for (celltype in unique(ranking$cell_type)) {
    one <- ranking[ranking$cell_type == celltype, , drop = FALSE]
    conditions <- unique(as.character(one$condition))
    if (length(conditions) < 2L) next
    condition_pairs <- utils::combn(conditions, 2L, simplify = FALSE)
    for (pair in condition_pairs) {
      a <- one[one$condition == pair[[1L]], , drop = FALSE]
      b <- one[one$condition == pair[[2L]], , drop = FALSE]
      b <- b[match(a$row_id, b$row_id), , drop = FALSE]
      if (anyNA(b$row_id)) {
        stop("Condition ranking tables contain different reaction targets.",
             call. = FALSE)
      }
      contrast_rows[[contrast_index]] <- data.frame(
        row_id = a$row_id,
        reaction_id = a$reaction_id,
        target_direction = a$target_direction,
        medium_scenario = a$medium_scenario,
        cell_type = celltype,
        condition_a = pair[[1L]],
        condition_b = pair[[2L]],
        median_penalty_a = a$median_penalty,
        median_penalty_b = b$median_penalty,
        priority_rank_a = a$priority_rank,
        priority_rank_b = b$priority_rank,
        delta_support_b_minus_a = b$support_score - a$support_score,
        higher_supported_condition = ifelse(
          b$support_score > a$support_score,
          pair[[2L]],
          ifelse(a$support_score > b$support_score, pair[[1L]], "tie")
        ),
        descriptive_only = TRUE,
        biological_replicate_inference = FALSE,
        stringsAsFactors = FALSE
      )
      contrast_index <- contrast_index + 1L
    }
  }
  contrast <- if (length(contrast_rows)) {
    do.call(rbind, contrast_rows)
  } else {
    data.frame()
  }
  if (nrow(contrast)) rownames(contrast) <- NULL

  n_conditions <- length(unique(as.character(ranking$condition)))
  analysis_mode <- if (n_conditions == 1L) {
    "single_condition_reaction_ranking"
  } else {
    "multi_condition_reaction_ranking_and_pairwise_comparison"
  }
  list(
    summary = ranking,
    ranking = ranking,
    contrast = contrast,
    analysis_mode = analysis_mode,
    inference_policy = paste(
      "condition-pooled metacells are descriptive pseudo-observations;",
      "reaction priority is ranked within each condition and cell type;",
      "biological-sample-level significance testing is not performed"
    )
  )
}
