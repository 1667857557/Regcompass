# Canonical Layer 1 implementation: condition-specific OOF edges enter the main path.

.rc_latent_metacell_expression <- function(
    counts, library_size, mu_min = 0.1,
    structural_zero_probability = 0.05, cell_type) {
  counts <- as.matrix(counts)
  library_size <- as.numeric(library_size)
  if (!is.numeric(counts) || is.null(rownames(counts)) ||
      is.null(colnames(counts)) || length(library_size) != ncol(counts) ||
      any(!is.finite(library_size)) || any(library_size <= 0) ||
      any(!is.finite(counts)) || any(counts < 0)) {
    stop("Latent RNA estimation requires aligned non-negative counts and depth.",
         call. = FALSE)
  }
  if (any(abs(counts - round(counts)) >
      sqrt(.Machine$double.eps) * pmax(1, abs(counts)))) {
    stop("Latent RNA estimation requires raw integer-like counts.", call. = FALSE)
  }
  if (!is.null(names(cell_type))) {
    cell_type <- as.character(cell_type[colnames(counts)])
  } else {
    cell_type <- as.character(cell_type)
  }
  if (length(cell_type) != ncol(counts) || anyNA(cell_type)) {
    stop("Cell-type labels do not align to metacell RNA counts.", call. = FALSE)
  }
  exposure <- library_size / 1e6
  observed_cpm <- sweep(counts, 2L, exposure, "/")
  types <- unique(cell_type)
  prior_shape <- prior_rate <- prior_mean <- prior_variance <- matrix(
    NA_real_, nrow(counts), length(types),
    dimnames = list(rownames(counts), types)
  )
  for (value in types) {
    x <- observed_cpm[, cell_type == value, drop = FALSE]
    mean_value <- rowMeans(x)
    variance_value <- if (ncol(x) > 1L) {
      apply(x, 1L, stats::var)
    } else {
      rep(NA_real_, nrow(x))
    }
    invalid <- !is.finite(variance_value) | variance_value <= 0
    variance_value[invalid] <- pmax(mean_value[invalid]^2 / 1e6, 1e-12)
    shape_value <- pmin(
      1e6, pmax(0.1, mean_value^2 / pmax(variance_value, 1e-12))
    )
    rate_value <- shape_value / pmax(mean_value, 1e-8)
    prior_mean[, value] <- mean_value
    prior_variance[, value] <- variance_value
    prior_shape[, value] <- shape_value
    prior_rate[, value] <- rate_value
  }
  shape_unit <- prior_shape[, cell_type, drop = FALSE]
  rate_unit <- prior_rate[, cell_type, drop = FALSE]
  colnames(shape_unit) <- colnames(rate_unit) <- colnames(counts)
  posterior_shape <- counts + shape_unit
  posterior_rate <- sweep(rate_unit, 2L, exposure, "+")
  latent_cpm <- posterior_shape / posterior_rate
  positive <- matrix(
    stats::pgamma(
      rep(mu_min, length(posterior_shape)),
      shape = as.numeric(posterior_shape),
      rate = as.numeric(posterior_rate),
      lower.tail = FALSE
    ),
    nrow = nrow(counts), dimnames = dimnames(counts)
  )
  zero_class <- matrix(
    "observed_positive", nrow(counts), ncol(counts), dimnames = dimnames(counts)
  )
  observed_zero <- counts == 0
  structural <- observed_zero & positive <= structural_zero_probability
  zero_class[observed_zero & !structural] <- "sampling_limited_zero"
  zero_class[structural] <- "credible_structural_zero"
  latent_cpm[structural] <- 0
  list(
    latent_cpm = latent_cpm,
    latent_log_expression = log1p(latent_cpm),
    posterior_positive_probability = positive,
    zero_class = zero_class,
    observed_zero = observed_zero,
    prior_mean = prior_mean,
    prior_variance = prior_variance,
    prior_shape = prior_shape,
    prior_rate = prior_rate,
    model = "gamma_poisson_empirical_bayes_by_cell_type_v2",
    prior_estimation_scope = "gene_by_cell_type"
  )
}

.rc_integrate_regulatory_support <- function(
    rna_support, regulatory_modifier, alpha = 1) {
  rna_support <- as.matrix(rna_support)
  regulatory_modifier <- as.matrix(regulatory_modifier)
  if (!identical(dimnames(rna_support), dimnames(regulatory_modifier))) {
    stop("RNA support and regulatory modifier matrices must align exactly.",
         call. = FALSE)
  }
  if (!isTRUE(all.equal(as.numeric(alpha), 1))) {
    stop("Canonical RegCompass requires `regulatory_alpha = 1`.", call. = FALSE)
  }
  C <- pmin(pmax(rna_support, 0), 1)
  fallback <- !is.finite(regulatory_modifier)
  R <- regulatory_modifier
  R[fallback] <- 0
  R <- pmin(pmax(R, -1), 1)
  multiplier <- 2^R
  out <- C * multiplier / (1 - C + C * multiplier)
  out[C <= 0] <- 0
  out[C >= 1] <- 1
  out[!is.finite(C)] <- NA_real_
  dimnames(out) <- dimnames(C)
  attr(out, "rna_only_fallback_mask") <- fallback
  attr(out, "integration_formula") <-
    "C_multiome=C_RNA*2^R/(1-C_RNA+C_RNA*2^R); nonfinite R:=0"
  out
}

.rc_projection_scale <- function(projection, unit_meta, celltype_col) {
  scale <- matrix(
    NA_real_, nrow(projection), ncol(projection), dimnames = dimnames(projection)
  )
  diagnostics <- list()
  for (gene in rownames(projection)) {
    for (celltype in unique(as.character(unit_meta[[celltype_col]]))) {
      units <- unit_meta$unit_id[
        as.character(unit_meta[[celltype_col]]) == celltype
      ]
      value <- projection[gene, units]
      value <- value[is.finite(value)]
      robust <- if (length(value) >= 2L && diff(range(value)) > 1e-12) {
        max(
          stats::IQR(value) / 1.349,
          stats::mad(value, constant = 1.4826),
          sqrt(mean(value^2)), 1e-6, na.rm = TRUE
        )
      } else {
        1
      }
      scale[gene, units] <- robust
      diagnostics[[length(diagnostics) + 1L]] <- data.frame(
        target = gene,
        cell_type = celltype,
        projection_scale = robust,
        stringsAsFactors = FALSE
      )
    }
  }
  list(scale = scale, diagnostics = do.call(rbind, diagnostics))
}

.rc_scaled_oof_modifier <- function(projection, reliability, scale) {
  if (!identical(dimnames(projection), dimnames(reliability)) ||
      !identical(dimnames(projection), dimnames(scale))) {
    stop("Projection, reliability, and calibration scale must align.",
         call. = FALSE)
  }
  value <- reliability * tanh(projection / scale)
  value[!is.finite(projection) | !is.finite(reliability)] <- NA_real_
  dimnames(value) <- dimnames(projection)
  value
}

.rc_gpr_best_group_fraction <- function(parsed_gpr, gene_available) {
  gene_available <- as.matrix(gene_available)
  answer <- matrix(
    NA_real_, length(parsed_gpr), ncol(gene_available),
    dimnames = list(names(parsed_gpr), colnames(gene_available))
  )
  genes_available <- tolower(rownames(gene_available))
  for (reaction_id in names(parsed_gpr)) {
    groups <- parsed_gpr[[reaction_id]]
    if (!length(groups)) next
    fractions <- vapply(groups, function(group) {
      requested <- unique(tolower(group))
      index <- match(requested, genes_available)
      vapply(seq_len(ncol(gene_available)), function(unit) {
        observed <- rep(FALSE, length(requested))
        mapped <- !is.na(index)
        observed[mapped] <- gene_available[index[mapped], unit]
        mean(observed)
      }, numeric(1))
    }, numeric(ncol(gene_available)))
    if (is.null(dim(fractions))) {
      fractions <- matrix(fractions, nrow = ncol(gene_available))
    }
    answer[reaction_id, ] <- apply(fractions, 1L, max)
  }
  answer
}

.rc_condition_pando_primary_projection <- function(
    grn_result, membership, unit_meta, genes) {
  primary <- reliability <- matrix(
    NA_real_, length(genes), nrow(unit_meta),
    dimnames = list(genes, unit_meta$unit_id)
  )
  coverage <- list()
  for (fit in grn_result$condition_grn_fits) {
    .rc_require_pando_condition_grn_fit(fit)
    cell <- Pando::project_condition_grn_primary_cells(
      object = grn_result$pando_grn_data,
      fit = fit,
      scale = "std",
      targets = intersect(genes, tolower(unique(fit$edge_table$target))),
      nonestimable = "structural_zero"
    )
    projected <- Pando::aggregate_condition_grn_projection(
      cell, membership, group_col = "metacell_id"
    )
    score <- t(as.matrix(projected$gene_score))
    rownames(score) <- tolower(rownames(score))
    target <- intersect(rownames(score), rownames(primary))
    units <- intersect(colnames(score), colnames(primary))
    primary[target, units] <- score[target, units, drop = FALSE]
    pooled <- fit$target_rsq_oof_pooled
    names(pooled) <- tolower(names(pooled))
    target <- intersect(names(pooled), genes)
    q <- sqrt(pmin(1, pmax(0, as.numeric(pooled[target]))))
    available <- fit$predictive_oof_available[
      match(target, tolower(names(fit$predictive_oof_available)))
    ] & fit$oof_cell_coverage[
      match(target, tolower(names(fit$oof_cell_coverage)))
    ] == 1
    q[!available] <- NA_real_
    units <- rownames(projected$gene_score)
    reliability[target, units] <- matrix(
      q, nrow = length(target), ncol = length(units)
    )
    coverage[[length(coverage) + 1L]] <-
      projected$source_projection$target_condition_status
  }
  list(
    primary = primary,
    reliability = reliability,
    coverage = .rc_bind_frames_fill(coverage),
    origin = "outer_condition_stratified_cell_oof",
    pando_schema = "pando_condition_grn_primary_projection_v1",
    condition_unique_edges_allowed = TRUE,
    common_support_required = FALSE
  )
}

.rc_cell_first_projection_layer1 <- function(
    grn_result, metacell_object, membership, metacell_meta, gem,
    condition_col, celltype_col, rna_assay,
    projection_component = "condition", regulatory_alpha = 1,
    gpr_and_method = "min", gene_half_saturation = 1,
    parallel = TRUE, BPPARAM = NULL) {
  if (!identical(projection_component, "condition")) {
    stop("Only the regulatory projection can enter Layer 1.", call. = FALSE)
  }
  if (!isTRUE(all.equal(as.numeric(regulatory_alpha), 1))) {
    stop("Canonical RegCompass requires `regulatory_alpha = 1`.", call. = FALSE)
  }
  if (!is.data.frame(membership) ||
      !all(c("cell_id", "metacell_id") %in% colnames(membership)) ||
      anyDuplicated(membership$cell_id)) {
    stop("SuperCell membership must map every cell exactly once.", call. = FALSE)
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
    stop("Metacell metadata do not align to the metacell object.", call. = FALSE)
  }
  parsed <- rc_parse_gpr_table(gem$gpr_table)
  gpr_genes <- unique(tolower(unlist(parsed, use.names = FALSE)))
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
  mode <- grn_result$analysis_mode %||% "condition_grn"
  projection <- if (identical(mode, "standard_pando")) {
    standard <- .rc_standard_pando_projection(
      grn_result, membership, unit_meta, condition_col, celltype_col,
      rna_assay = grn_result$rna_assay %||% "RNA",
      atac_assay = grn_result$atac_assay %||% "ATAC",
      target_genes = genes
    )
    list(
      primary = standard$projection,
      reliability = standard$reliability,
      coverage = standard$coverage,
      origin = standard$projection_origin,
      pando_schema = "standard_pando_network",
      condition_unique_edges_allowed = FALSE,
      common_support_required = FALSE
    )
  } else {
    .rc_condition_pando_primary_projection(
      grn_result, membership, unit_meta, genes
    )
  }
  calibration <- .rc_projection_scale(
    projection$primary, unit_meta, celltype_col
  )
  modifier <- .rc_scaled_oof_modifier(
    projection$primary, projection$reliability, calibration$scale
  )
  gene_support_rna <- rc_gene_score(
    latent$latent_log_expression,
    mode = "absolute",
    half_saturation = gene_half_saturation
  )
  gene_support_multiome <- .rc_integrate_regulatory_support(
    gene_support_rna, modifier, alpha = 1
  )
  reaction_rna <- rc_reaction_capacity(
    parsed, gene_support_rna,
    promiscuity_mode = "none",
    and_method = gpr_and_method,
    or_method = "sum",
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )
  reaction_primary <- rc_reaction_capacity(
    parsed, gene_support_multiome,
    promiscuity_mode = "none",
    and_method = gpr_and_method,
    or_method = "sum",
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )
  support_fraction <- .rc_gpr_best_group_fraction(
    parsed, is.finite(modifier)
  )
  fallback <- !is.finite(modifier)
  list(
    schema_version = "regcompass_regulatory_layer1_v2",
    analysis_mode = mode,
    reaction_expression = reaction_primary,
    reaction_expression_condition_oof = reaction_primary,
    reaction_expression_rna_only = reaction_rna,
    reaction_expression_available = is.finite(reaction_primary),
    rna_metacell_latent_log_expression = latent$latent_log_expression,
    rna_metacell_latent_cpm = latent$latent_cpm,
    posterior_positive_probability = latent$posterior_positive_probability,
    posterior_zero_probability = 1 - latent$posterior_positive_probability,
    rna_zero_class = latent$zero_class,
    gene_support_rna = gene_support_rna,
    gene_support_multiome = gene_support_multiome,
    gene_projection_condition_oof = projection$primary,
    gene_projection_raw = projection$primary,
    gene_projection_scale = calibration$scale,
    gene_regulatory_reliability = projection$reliability,
    gene_regulatory_reliability_available = is.finite(projection$reliability),
    gene_regulatory_available = is.finite(projection$primary),
    gene_regulatory_modifier = modifier,
    projection_coverage = projection$coverage,
    projection_calibration = calibration$diagnostics,
    reaction_condition_support_fraction = support_fraction,
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
      prior_estimation_scope = latent$prior_estimation_scope
    ),
    zero_diagnostics = list(
      observed_zero_fraction = rowMeans(latent$observed_zero),
      posterior_sampling_zero_fraction = rowMeans(
        latent$zero_class == "sampling_limited_zero"
      ),
      credible_structural_zero_fraction = rowMeans(
        latent$zero_class == "credible_structural_zero"
      )
    ),
    capacity_params = list(
      regulatory_alpha = 1,
      regulatory_odds_budget = 2,
      gene_half_saturation = gene_half_saturation,
      regulatory_mode = mode,
      link_function = "tanh(G/shared_scale)",
      promiscuity_mode = "none",
      and_method = gpr_and_method,
      or_method = "sum",
      parallel = parallel
    ),
    projection_provenance = list(
      analysis_mode = mode,
      pando_schema = projection$pando_schema,
      projection_origin = projection$origin,
      projection_used_for_penalty = TRUE,
      full_fit_projection_used_for_penalty =
        identical(mode, "standard_pando"),
      condition_coefficients_calculated = identical(mode, "condition_grn"),
      primary_support_policy = if (identical(mode, "condition_grn")) {
        "condition_estimable"
      } else {
        "standard_pando_full_fit"
      },
      condition_unique_edges_allowed =
        isTRUE(projection$condition_unique_edges_allowed),
      common_support_required = FALSE,
      supercell_membership = "membership_table(cell_id, metacell_id)",
      unavailable_target_policy = "rna_only_neutral_modifier_fallback",
      nonestimable_edge_policy = if (identical(mode, "condition_grn")) {
        "condition_specific_structural_zero_enters_main_analysis"
      } else {
        "not_applicable_standard_pando"
      }
    ),
    inference_class = "metacell_statistical_unit_within_dataset",
    statistical_unit = "metacell",
    metacell_statistical_inference = TRUE,
    biological_replicate_inference = FALSE
  )
}
