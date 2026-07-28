.rc_build_condition_pooled_layer1 <- function(
    metacell_object, meta_modules, gem, metacell_meta,
    condition_col = "condition", celltype_col = "cell_type",
    rna_assay = "RNA", atac_assay = "ATAC",
    regulatory_alpha = 1,
    gpr_and_method = c("min", "median", "mean"),
    gene_half_saturation = getOption("RegCompassR.cpm_half_saturation", 1),
    parallel = TRUE, BPPARAM = NULL) {
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("`parallel` must be TRUE or FALSE.", call. = FALSE)
  }
  gpr_and_method <- match.arg(gpr_and_method)
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
    stop("Duplicated GPR gene identifiers after case normalization.",
         call. = FALSE)
  }

  unit_meta <- metacell_meta
  id_col <- if ("metacell_id" %in% colnames(unit_meta)) {
    "metacell_id"
  } else if ("pool_id" %in% colnames(unit_meta)) {
    "pool_id"
  } else {
    stop("Metacell metadata lack metacell_id/pool_id.", call. = FALSE)
  }
  required_meta <- c(condition_col, celltype_col)
  missing_meta <- setdiff(required_meta, colnames(unit_meta))
  if (length(missing_meta)) {
    stop(
      "Metacell metadata lack condition/cell-type columns: ",
      paste(missing_meta, collapse = ", "), call. = FALSE
    )
  }
  if ("n_celltypes" %in% colnames(unit_meta) &&
      any(unit_meta$n_celltypes != 1L, na.rm = TRUE)) {
    stop("Layer 1 requires SuperCell2 label-pure metacells.", call. = FALSE)
  }
  if ("mixed_celltype_metacell" %in% colnames(unit_meta) &&
      any(unit_meta$mixed_celltype_metacell %in% TRUE, na.rm = TRUE)) {
    stop("Layer 1 cannot use mixed-cell-type metacells.", call. = FALSE)
  }
  unit_meta$pool_id <- as.character(unit_meta[[id_col]])
  unit_meta$unit_id <- unit_meta$pool_id
  unit_meta <- unit_meta[
    match(colnames(rna_logcpm), unit_meta$pool_id), , drop = FALSE
  ]
  if (anyNA(unit_meta$pool_id)) {
    stop("Metacell metadata do not align with RNA counts.", call. = FALSE)
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
  gene_multiome_support <- .rc_integrate_regulatory_support(
    gene_rna_support,
    modifier,
    alpha = regulatory_alpha
  )
  reaction_expression <- rc_reaction_capacity(
    parsed,
    gene_multiome_support,
    promiscuity_mode = "none",
    and_method = gpr_and_method,
    or_method = "sum",
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )

  list(
    schema_version = "regcompass_condition_only_layer1_v5_supercell_label",
    reaction_expression = reaction_expression,
    rna_metacell_logcpm = rna_logcpm,
    gene_support_rna = gene_rna_support,
    gene_regulatory_modifier = modifier,
    gene_support_multiome = gene_multiome_support,
    parsed_gpr = parsed,
    gpr_diagnostics = rc_gpr_diagnostics(parsed, rownames(rna_logcpm)),
    unit_meta = unit_meta,
    metacell_meta = unit_meta,
    layer1_unit = "condition_stratified_supercell2_label_exact_celltype",
    capacity_params = list(
      regulatory_alpha = regulatory_alpha,
      gene_half_saturation = gene_half_saturation,
      promiscuity_mode = "none",
      and_method = gpr_and_method,
      or_method = "sum",
      parallel = parallel,
      bpparam_class = if (is.null(BPPARAM)) {
        "auto_or_sequential"
      } else if (identical(BPPARAM, FALSE)) {
        "sequential"
      } else {
        class(BPPARAM)[[1L]]
      }
    ),
    evidence_formula = paste(
      "bootstrap-stable multitask coefficient-weighted ATAC modifier ->",
      "zero-preserving RNA support log-odds update ->",
      paste0("COMPASS-compatible ", gpr_and_method, " GPR-AND"),
      "and additive isozyme OR"
    ),
    evidence_inputs = c(
      "target_gene_RNA",
      "peak_ATAC",
      "condition_x_celltype_bootstrap_stable_coefficients"
    )
  )
}
