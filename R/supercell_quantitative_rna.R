.rc_pando_rna_objects <- function(grn_result) {
  condition_objects <- grn_result$pando_grn_data_by_cell_type %||% list()
  if (!length(condition_objects) &&
      inherits(grn_result$pando_grn_data, "GRNData")) {
    condition_objects <- list(condition_legacy = grn_result$pando_grn_data)
  }
  standard_objects <- grn_result$standard_pando_objects %||% list()
  objects <- c(condition_objects, standard_objects)
  objects <- objects[!vapply(objects, is.null, logical(1))]
  if (!length(objects)) {
    stop(
      "Quantitative RNA aggregation requires the cell-level Pando objects retained by Stage 1.",
      call. = FALSE
    )
  }
  if (is.null(names(objects)) || anyNA(names(objects)) ||
      any(!nzchar(names(objects))) || anyDuplicated(names(objects))) {
    stop("Pando RNA source object names must be unique and non-empty.",
         call. = FALSE)
  }
  invalid <- !vapply(objects, inherits, logical(1), what = "GRNData")
  if (any(invalid)) {
    stop(
      "Quantitative RNA aggregation requires GRNData objects for every routed cell type.",
      call. = FALSE
    )
  }
  objects
}

.rc_validate_pando_rna_cell_partition <- function(
    source_cells, membership_cells) {
  if (!is.list(source_cells) || !length(source_cells) ||
      is.null(names(source_cells)) || anyNA(names(source_cells)) ||
      any(!nzchar(names(source_cells))) || anyDuplicated(names(source_cells))) {
    stop("Pando RNA source cell sets require unique non-empty source names.",
         call. = FALSE)
  }
  source_cells <- lapply(source_cells, function(cells) {
    cells <- as.character(cells)
    if (!length(cells) || anyNA(cells) || any(!nzchar(cells)) ||
        anyDuplicated(cells)) {
      stop("Every Pando RNA source requires unique non-empty cell IDs.",
           call. = FALSE)
    }
    cells
  })
  membership_cells <- as.character(membership_cells)
  if (!length(membership_cells) || anyNA(membership_cells) ||
      any(!nzchar(membership_cells)) || anyDuplicated(membership_cells)) {
    stop("SuperCell membership requires unique non-empty cell IDs.",
         call. = FALSE)
  }
  observed <- unlist(source_cells, use.names = FALSE)
  if (anyDuplicated(observed)) {
    duplicated_cells <- unique(observed[duplicated(observed)])
    stop(
      "A cell occurs in more than one routed Pando RNA source; first duplicated cell: ",
      duplicated_cells[[1L]], ".", call. = FALSE
    )
  }
  missing <- setdiff(membership_cells, observed)
  extra <- setdiff(observed, membership_cells)
  if (length(missing) || length(extra)) {
    stop(
      "Stage 1 Pando RNA sources and SuperCell membership are not an exact cell partition; ",
      "missing=", length(missing), ", extra=", length(extra), ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.rc_single_cell_linear_cpm <- function(
    counts, genes, scale_factor = 1e6) {
  if (is.null(dim(counts)) || is.null(rownames(counts)) ||
      is.null(colnames(counts)) || anyDuplicated(rownames(counts)) ||
      anyDuplicated(colnames(counts))) {
    stop("Single-cell RNA counts require unique gene and cell IDs.",
         call. = FALSE)
  }
  if (!is.numeric(scale_factor) || length(scale_factor) != 1L ||
      !is.finite(scale_factor) || scale_factor <= 0) {
    stop("`scale_factor` must be one positive finite number.", call. = FALSE)
  }
  genes <- unique(tolower(trimws(as.character(genes))))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  feature_key <- tolower(rownames(counts))
  if (anyDuplicated(feature_key)) {
    stop("Duplicated RNA genes after case normalization.", call. = FALSE)
  }
  index <- match(genes, feature_key)
  if (anyNA(index)) {
    missing <- genes[is.na(index)]
    stop(
      "Cell-level Pando RNA counts are missing GPR gene(s); first missing gene: ",
      missing[[1L]], ".",
      call. = FALSE
    )
  }
  library_size <- as.numeric(Matrix::colSums(counts))
  if (any(!is.finite(library_size)) || any(library_size <= 0)) {
    stop("Every cell must have a positive finite RNA library size.",
         call. = FALSE)
  }
  selected <- counts[index, , drop = FALSE]
  rownames(selected) <- genes
  expression <- selected %*% Matrix::Diagonal(
    x = scale_factor / library_size
  )
  dimnames(expression) <- list(genes, colnames(counts))
  list(
    expression = expression,
    library_size = stats::setNames(library_size, colnames(counts)),
    scale_factor = scale_factor
  )
}

.rc_equal_mean_supercell_expression <- function(
    cell_expression, membership, units) {
  if (!is.data.frame(membership) ||
      !all(c("cell_id", "metacell_id") %in% colnames(membership)) ||
      anyNA(membership$cell_id) || anyNA(membership$metacell_id) ||
      any(!nzchar(as.character(membership$cell_id))) ||
      any(!nzchar(as.character(membership$metacell_id))) ||
      anyDuplicated(membership$cell_id)) {
    stop("SuperCell membership must map every cell exactly once.",
         call. = FALSE)
  }
  if (is.null(dim(cell_expression)) || is.null(rownames(cell_expression)) ||
      is.null(colnames(cell_expression)) ||
      anyDuplicated(rownames(cell_expression)) ||
      anyDuplicated(colnames(cell_expression))) {
    stop("Cell-level quantitative RNA requires unique gene and cell IDs.",
         call. = FALSE)
  }
  units <- as.character(units)
  if (anyNA(units) || any(!nzchar(units)) || anyDuplicated(units)) {
    stop("Metacell unit IDs must be unique and non-empty.", call. = FALSE)
  }
  cell_ids <- as.character(membership$cell_id)
  if (!setequal(colnames(cell_expression), cell_ids)) {
    stop(
      "Cell-level quantitative RNA and SuperCell membership must cover exactly the same cells.",
      call. = FALSE
    )
  }
  membership <- membership[
    match(colnames(cell_expression), cell_ids), , drop = FALSE
  ]
  unit_index <- match(as.character(membership$metacell_id), units)
  if (anyNA(unit_index)) {
    stop("SuperCell membership contains metacell IDs absent from Layer 1 units.",
         call. = FALSE)
  }
  n_cells <- tabulate(unit_index, nbins = length(units))
  if (any(n_cells < 1L)) {
    stop("Every Layer 1 metacell must contain at least one member cell.",
         call. = FALSE)
  }
  averaging <- Matrix::sparseMatrix(
    i = seq_len(nrow(membership)),
    j = unit_index,
    x = 1 / n_cells[unit_index],
    dims = c(nrow(membership), length(units)),
    dimnames = list(colnames(cell_expression), units)
  )
  expression <- as.matrix(cell_expression %*% averaging)
  dimnames(expression) <- list(rownames(cell_expression), units)
  list(
    expression = expression,
    n_cells = stats::setNames(as.integer(n_cells), units),
    aggregation = "equal_mean_over_single_cell_linear_cpm",
    library_size_weighted = FALSE
  )
}

.rc_quantitative_supercell_rna <- function(
    grn_result, membership, units, genes, rna_assay,
    scale_factor = 1e6) {
  objects <- .rc_pando_rna_objects(grn_result)
  membership_cells <- as.character(membership$cell_id)
  object_cells <- lapply(objects, function(object) {
    colnames(object@data)
  })
  .rc_validate_pando_rna_cell_partition(object_cells, membership_cells)
  expression_parts <- list()
  library_parts <- list()
  source_rows <- list()
  for (name in names(objects)) {
    object <- objects[[name]]
    cells <- colnames(object@data)
    counts <- .rc_get_assay_counts(object@data, rna_assay)
    if (!all(cells %in% colnames(counts))) {
      stop("Pando RNA count columns do not match the stored cell IDs.",
           call. = FALSE)
    }
    normalized <- .rc_single_cell_linear_cpm(
      counts[, cells, drop = FALSE], genes, scale_factor = scale_factor
    )
    expression_parts[[name]] <- normalized$expression
    library_parts[[name]] <- normalized$library_size
    source_rows[[name]] <- data.frame(
      source = name,
      n_cells = length(cells),
      stringsAsFactors = FALSE
    )
  }
  if (!length(expression_parts)) {
    stop("No Stage 1 Pando cells overlap the SuperCell membership.",
         call. = FALSE)
  }
  observed_cells <- unlist(lapply(expression_parts, colnames), use.names = FALSE)
  if (anyDuplicated(observed_cells)) {
    duplicated_cells <- unique(observed_cells[duplicated(observed_cells)])
    stop(
      "A cell occurs in more than one routed Pando RNA source; first duplicated cell: ",
      duplicated_cells[[1L]], ".",
      call. = FALSE
    )
  }
  cell_expression <- do.call(cbind, expression_parts)
  cell_expression <- cell_expression[, membership_cells, drop = FALSE]
  averaged <- .rc_equal_mean_supercell_expression(
    cell_expression, membership, units
  )
  library_size <- unlist(library_parts, use.names = FALSE)
  names(library_size) <- observed_cells
  library_size <- library_size[membership_cells]
  if (anyNA(library_size)) {
    stop("Cell-level RNA library-size provenance is incomplete.", call. = FALSE)
  }
  list(
    expression = averaged$expression,
    n_cells = averaged$n_cells,
    cell_library_size = library_size,
    scale_factor = scale_factor,
    aggregation = averaged$aggregation,
    library_size_weighted = averaged$library_size_weighted,
    source_summary = do.call(rbind, source_rows),
    model = "supercell_equal_mean_single_cell_linear_cpm_v1"
  )
}

.rc_cell_first_projection_layer1_v6 <- function(
    grn_result, metacell_object, membership, metacell_meta, gem,
    condition_col, celltype_col, rna_assay,
    gpr_and_method = "min", gene_half_saturation = 1,
    parallel = TRUE, BPPARAM = NULL) {
  if (!is.data.frame(membership) ||
      !all(c("cell_id", "metacell_id") %in% colnames(membership)) ||
      anyDuplicated(membership$cell_id)) {
    stop("SuperCell membership must map every cell exactly once.",
         call. = FALSE)
  }
  id_col <- if ("metacell_id" %in% colnames(metacell_meta)) {
    "metacell_id"
  } else {
    "pool_id"
  }
  unit_meta <- as.data.frame(metacell_meta)
  unit_meta$unit_id <- as.character(unit_meta[[id_col]])
  unit_meta$pool_id <- unit_meta$unit_id
  units <- colnames(metacell_object)
  unit_meta <- unit_meta[match(units, unit_meta$unit_id), , drop = FALSE]
  if (anyNA(unit_meta$unit_id)) {
    stop("Metacell metadata do not align to the metacell object.",
         call. = FALSE)
  }
  parsed <- rc_parse_gpr_table(gem$gpr_table)
  gpr_genes <- unique(tolower(unlist(parsed, use.names = FALSE)))

  # Aggregated raw metacell counts remain the input to the legacy bounded
  # structural-confidence model only. They are no longer the quantitative LP
  # expression estimator.
  counts <- .rc_get_assay_counts(metacell_object, rna_assay)
  library_size <- Matrix::colSums(counts)
  rna_counts <- counts[
    tolower(rownames(counts)) %in% gpr_genes, units, drop = FALSE
  ]
  rownames(rna_counts) <- tolower(rownames(rna_counts))
  if (anyDuplicated(rownames(rna_counts))) {
    stop("Duplicated GPR genes after case normalization.", call. = FALSE)
  }
  cell_type <- stats::setNames(
    as.character(unit_meta[[celltype_col]]), unit_meta$unit_id
  )
  latent <- .rc_latent_metacell_expression(
    rna_counts, library_size[units], cell_type = cell_type
  )
  genes <- rownames(latent$latent_log_expression)

  # Quantitative COMPASS path follows the SuperCell representative-state
  # definition. Each original cell is normalized on the linear CPM scale using
  # its own complete RNA library, then member cells are averaged with equal
  # weight inside the exact final SuperCell membership. This is deliberately not
  # CPM(sum(counts)), which would weight cells by library size.
  quantitative_rna <- .rc_quantitative_supercell_rna(
    grn_result = grn_result,
    membership = membership,
    units = units,
    genes = genes,
    rna_assay = rna_assay,
    scale_factor = 1e6
  )
  if (!identical(
        dimnames(quantitative_rna$expression),
        list(genes, units)
      )) {
    stop("Quantitative SuperCell RNA is misaligned to Layer 1 genes or units.",
         call. = FALSE)
  }

  mode <- grn_result$analysis_mode %||% "condition_grn"
  projection <- .rc_project_pando_by_celltype(
    grn_result = grn_result,
    membership = membership,
    unit_meta = unit_meta,
    genes = genes,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = grn_result$rna_assay %||% "RNA",
    atac_assay = grn_result$atac_assay %||% "ATAC"
  )
  calibration <- .rc_projection_scale(
    projection$projection, unit_meta, celltype_col
  )
  modifier <- .rc_scaled_regulatory_modifier(
    projection$projection, projection$reliability, calibration$scale
  )

  gene_expression_quantitative_rna <- quantitative_rna$expression
  gene_expression_quantitative_multiome <-
    .rc_integrate_regulatory_expression(
      gene_expression_quantitative_rna, modifier
    )
  attr(
    gene_expression_quantitative_multiome,
    "integration_formula"
  ) <- paste(
    "X_multiome=X_RNA*2^R;",
    "X_RNA=equal_mean_single_cell_linear_CPM_by_SuperCell_membership;",
    "nonfinite R:=0; R clipped to [-1,1]"
  )
  reaction_quantitative_rna <- rc_reaction_capacity(
    parsed, gene_expression_quantitative_rna,
    promiscuity_mode = "none",
    and_method = gpr_and_method, or_method = "sum",
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )
  reaction_quantitative_multiome <- rc_reaction_capacity(
    parsed, gene_expression_quantitative_multiome,
    promiscuity_mode = "none",
    and_method = gpr_and_method, or_method = "sum",
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )

  # Bounded support intentionally keeps the previous latent-CPM model for
  # CORDA2/structural confidence in this change. The quantitative LP route above
  # is completely independent of that empirical-Bayes shrinkage.
  gene_support_rna <- rc_gene_score(
    latent$latent_log_expression,
    mode = "absolute",
    half_saturation = gene_half_saturation
  )
  gene_support_multiome <- .rc_integrate_regulatory_support(
    gene_support_rna, modifier, alpha = 1
  )
  reaction_structural_rna <- rc_reaction_capacity(
    parsed, gene_support_rna, promiscuity_mode = "none",
    and_method = gpr_and_method, or_method = "sum",
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )
  reaction_structural_multiome <- rc_reaction_capacity(
    parsed, gene_support_multiome, promiscuity_mode = "none",
    and_method = gpr_and_method, or_method = "sum",
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )
  attr(
    reaction_structural_multiome,
    "regcompass_quantitative_penalty_route"
  ) <- "multiome"
  attr(
    reaction_structural_rna,
    "regcompass_quantitative_penalty_route"
  ) <- "rna_only"

  support_fraction <- .rc_gpr_best_group_fraction(
    parsed, is.finite(modifier)
  )
  fallback <- !is.finite(modifier)
  cell_library <- as.numeric(quantitative_rna$cell_library_size)
  list(
    schema_version = "regcompass_regulatory_layer1_v6",
    analysis_mode = mode,

    reaction_expression = reaction_structural_multiome,
    reaction_expression_rna_only = reaction_structural_rna,
    reaction_expression_available = is.finite(reaction_structural_multiome),
    reaction_structural_support = reaction_structural_multiome,
    reaction_structural_support_rna_only = reaction_structural_rna,

    reaction_expression_quantitative = reaction_quantitative_multiome,
    reaction_expression_quantitative_rna_only = reaction_quantitative_rna,
    reaction_expression_quantitative_available =
      is.finite(reaction_quantitative_multiome),

    rna_metacell_mean_single_cell_cpm = gene_expression_quantitative_rna,
    rna_metacell_latent_log_expression = latent$latent_log_expression,
    rna_metacell_latent_cpm = latent$latent_cpm,
    posterior_positive_probability = latent$posterior_positive_probability,
    posterior_zero_probability = latent$posterior_zero_probability,
    rna_zero_class = latent$zero_class,
    eb_prior_weight = latent$prior_weight,
    eb_observation_weight = latent$observation_weight,
    gene_expression_quantitative_rna = gene_expression_quantitative_rna,
    gene_expression_quantitative_multiome =
      gene_expression_quantitative_multiome,
    gene_support_rna = gene_support_rna,
    gene_support_multiome = gene_support_multiome,
    gene_projection = projection$projection,
    gene_projection_scale = calibration$scale,
    gene_regulatory_reliability = projection$reliability,
    gene_regulatory_reliability_available =
      is.finite(projection$reliability),
    gene_regulatory_available = is.finite(projection$projection),
    gene_regulatory_modifier = modifier,
    projection_coverage = projection$coverage,
    projection_calibration = calibration$diagnostics,
    reaction_regulatory_support_fraction = support_fraction,
    parsed_gpr = parsed,
    gpr_diagnostics = rc_gpr_diagnostics(parsed, genes),
    unit_meta = unit_meta,
    metacell_meta = unit_meta,
    layer1_unit = "native_SuperCell_metacell",
    regulatory_fallback = list(
      policy = "rna_only_for_nonfinite_pando_modifier",
      neutral_modifier = 0,
      gene_metacell_mask = fallback,
      n_fallback = sum(fallback),
      fallback_fraction = mean(fallback)
    ),
    depth_diagnostics = list(
      rna_library_size = stats::setNames(
        as.numeric(library_size[units]), units
      ),
      latent_expression_model = latent$model,
      prior_estimation_scope = latent$prior_estimation_scope,
      posterior_update_scope = latent$posterior_update_scope,
      eb_weight_library_size_dependent = FALSE,
      quantitative_expression_model = quantitative_rna$model,
      quantitative_normalization_scale = quantitative_rna$scale_factor,
      quantitative_aggregation = quantitative_rna$aggregation,
      quantitative_library_size_weighted =
        quantitative_rna$library_size_weighted,
      quantitative_metacell_n_cells = quantitative_rna$n_cells,
      quantitative_cell_library_size_summary = c(
        min = min(cell_library),
        median = stats::median(cell_library),
        max = max(cell_library)
      ),
      quantitative_rna_source_summary = quantitative_rna$source_summary
    ),
    zero_diagnostics = list(
      observed_zero_fraction = rowMeans(latent$observed_zero),
      mean_posterior_zero_probability = rowMeans(
        latent$posterior_zero_probability
      ),
      maximum_posterior_zero_probability = apply(
        latent$posterior_zero_probability, 1L, max
      )
    ),
    capacity_params = list(
      regulatory_odds_budget = 2,
      regulatory_expression_multiplier_budget = 2,
      gene_half_saturation = gene_half_saturation,
      regulatory_mode = mode,
      link_function = "tanh(G/shared_scale)",
      quantitative_gene_expression = "equal_mean_single_cell_linear_cpm",
      quantitative_regulatory_formula =
        "X_multiome=X_RNA*2^R; nonfinite R:=0",
      structural_gene_support =
        "log1p(latent_cpm)/(log1p(latent_cpm)+gene_half_saturation)",
      structural_regulatory_formula =
        "C_multiome=C_RNA*2^R/(1-C_RNA+C_RNA*2^R)",
      promiscuity_mode = "none",
      and_method = gpr_and_method,
      or_method = "sum",
      parallel = parallel
    ),
    quantitative_penalty_contract = list(
      baseline_gene_expression = "equal_mean_single_cell_linear_cpm",
      single_cell_normalization =
        "raw_counts_per_cell_divided_by_complete_cell_RNA_library_times_1e6",
      metacell_aggregation = "equal_mean_by_exact_SuperCell_membership",
      library_size_weighted_metacell_average = FALSE,
      regulatory_multiplier = "2^R",
      regulatory_modifier_range = c(-1, 1),
      gpr_and_method = gpr_and_method,
      gpr_or_method = "sum",
      penalty_formula = "1/(1+log2(1+max(E_quantitative,0)))",
      structural_route_marker = "regcompass_quantitative_penalty_route",
      primary_route = "multiome",
      rna_control_route = "rna_only",
      bounded_support_excluded_from_lp_penalty = TRUE
    ),
    structural_support_contract = list(
      gene_support_range = c(0, 1),
      gene_half_saturation = gene_half_saturation,
      regulatory_update = "bounded_odds",
      intended_use = "CORDA2_and_structural_confidence",
      quantitative_lp_penalty = FALSE,
      latent_cpm_structural_only = TRUE
    ),
    projection_provenance = list(
      analysis_mode = mode,
      pando_schema = projection$pando_schema,
      projection_origin = projection$origin,
      projection_used_for_penalty = TRUE,
      projection_name = projection$projection_name,
      condition_coefficients_calculated =
        isTRUE(projection$condition_coefficients_calculated),
      cell_type_analysis_mode = projection$cell_type_analysis_mode,
      supercell_membership = "membership_table(cell_id, metacell_id)",
      quantitative_rna_source = "Pando_retained_cell_level_raw_RNA_counts",
      quantitative_rna_aggregation =
        "equal_mean_after_per_cell_linear_CPM_normalization",
      unavailable_target_policy = "rna_only_neutral_modifier_fallback",
      nonestimable_edge_policy = projection$nonestimable_policy
    ),
    inference_class = "metacell_statistical_unit_within_dataset",
    statistical_unit = "metacell",
    metacell_statistical_inference = TRUE,
    biological_replicate_inference = FALSE
  )
}
