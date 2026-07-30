# Cell-type-specific RNA priors and audited RNA-only regulatory fallback.

.rc_layer1_policy_state <- new.env(parent = emptyenv())

.rc_latent_metacell_expression <- function(
    counts, library_size, mu_min = 0.1,
    structural_zero_probability = 0.05,
    cell_type = NULL) {
  counts <- as.matrix(counts)
  library_size <- as.numeric(library_size)
  if (!is.numeric(counts) || is.null(rownames(counts)) ||
      is.null(colnames(counts)) ||
      length(library_size) != ncol(counts) ||
      any(!is.finite(library_size)) || any(library_size <= 0) ||
      any(!is.finite(counts)) || any(counts < 0)) {
    stop("Latent RNA estimation requires aligned non-negative counts and depth.",
         call. = FALSE)
  }
  integer_tolerance <- sqrt(.Machine$double.eps) * pmax(1, abs(counts))
  if (any(abs(counts - round(counts)) > integer_tolerance)) {
    stop(
      "Gamma-Poisson latent RNA estimation requires raw integer-like counts.",
      call. = FALSE
    )
  }
  if (is.null(cell_type)) {
    cell_type <- .rc_layer1_policy_state$cell_type
  }
  if (is.null(cell_type)) {
    stop(
      "Cell-type labels are required for metacell prior estimation.",
      call. = FALSE
    )
  }
  if (!is.null(names(cell_type))) {
    cell_type <- as.character(cell_type[colnames(counts)])
  } else {
    cell_type <- as.character(cell_type)
  }
  if (length(cell_type) != ncol(counts) || anyNA(cell_type) ||
      any(!nzchar(trimws(cell_type)))) {
    stop("Cell-type labels do not align to metacell RNA counts.",
         call. = FALSE)
  }

  exposure <- library_size / 1e6
  observed_cpm <- sweep(counts, 2L, exposure, "/")
  cell_types <- unique(cell_type)
  prior_shape_by_type <- matrix(
    NA_real_, nrow(counts), length(cell_types),
    dimnames = list(rownames(counts), cell_types)
  )
  prior_rate_by_type <- prior_shape_by_type
  prior_mean_by_type <- prior_shape_by_type
  prior_variance_by_type <- prior_shape_by_type

  for (value in cell_types) {
    selected <- cell_type == value
    x <- observed_cpm[, selected, drop = FALSE]
    prior_mean <- rowMeans(x)
    prior_variance <- if (ncol(x) > 1L) {
      apply(x, 1L, stats::var)
    } else {
      rep(NA_real_, nrow(x))
    }
    invalid_variance <- !is.finite(prior_variance) | prior_variance <= 0
    prior_variance[invalid_variance] <- pmax(
      prior_mean[invalid_variance]^2 / 1e6, 1e-12
    )
    prior_shape <- pmin(
      1e6, pmax(0.1, prior_mean^2 / pmax(prior_variance, 1e-12))
    )
    prior_rate <- prior_shape / pmax(prior_mean, 1e-8)
    prior_mean_by_type[, value] <- prior_mean
    prior_variance_by_type[, value] <- prior_variance
    prior_shape_by_type[, value] <- prior_shape
    prior_rate_by_type[, value] <- prior_rate
  }

  prior_shape_unit <- prior_shape_by_type[, cell_type, drop = FALSE]
  prior_rate_unit <- prior_rate_by_type[, cell_type, drop = FALSE]
  colnames(prior_shape_unit) <- colnames(counts)
  colnames(prior_rate_unit) <- colnames(counts)
  posterior_shape <- counts + prior_shape_unit
  posterior_rate <- sweep(prior_rate_unit, 2L, exposure, "+")
  latent_cpm <- posterior_shape / posterior_rate
  posterior_positive <- matrix(
    stats::pgamma(
      rep(mu_min, length(posterior_shape)),
      shape = as.numeric(posterior_shape),
      rate = as.numeric(posterior_rate),
      lower.tail = FALSE
    ),
    nrow = nrow(counts),
    dimnames = dimnames(counts)
  )
  zero_class <- matrix(
    "observed_positive", nrow(counts), ncol(counts),
    dimnames = dimnames(counts)
  )
  observed_zero <- counts == 0
  credible_structural <- observed_zero &
    posterior_positive <= structural_zero_probability
  zero_class[observed_zero & !credible_structural] <-
    "sampling_limited_zero"
  zero_class[credible_structural] <- "credible_structural_zero"
  latent_cpm[credible_structural] <- 0
  dimnames(latent_cpm) <- dimnames(counts)

  list(
    latent_cpm = latent_cpm,
    latent_log_expression = log1p(latent_cpm),
    posterior_positive_probability = posterior_positive,
    zero_class = zero_class,
    observed_zero = observed_zero,
    dispersion = 1 / prior_shape_by_type,
    prior_mean = prior_mean_by_type,
    prior_variance = prior_variance_by_type,
    prior_shape = prior_shape_by_type,
    prior_rate = prior_rate_by_type,
    prior_cell_type = stats::setNames(cell_type, colnames(counts)),
    prior_estimation_scope = "gene_by_cell_type",
    model = "gamma_poisson_empirical_bayes_by_cell_type_v2",
    marginal_count_model = "negative_binomial",
    mu_min = mu_min,
    structural_zero_probability = structural_zero_probability
  )
}

.rc_integrate_regulatory_support <- function(
    rna_support, regulatory_modifier, alpha = 1) {
  rna_support <- as.matrix(rna_support)
  regulatory_modifier <- as.matrix(regulatory_modifier)
  if (!identical(dim(rna_support), dim(regulatory_modifier)) ||
      !identical(dimnames(rna_support), dimnames(regulatory_modifier))) {
    stop("RNA support and regulatory modifier matrices must align exactly.",
         call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || !isTRUE(all.equal(as.numeric(alpha), 1))) {
    stop("Canonical RegCompass requires `regulatory_alpha = 1`.",
         call. = FALSE)
  }
  C <- pmin(pmax(rna_support, 0), 1)
  fallback <- !is.finite(regulatory_modifier)
  R <- regulatory_modifier
  R[fallback] <- 0
  R <- pmin(pmax(R, -1), 1)
  multiplier <- 2^R
  numerator <- C * multiplier
  denominator <- 1 - C + numerator
  out <- numerator / denominator
  out[C <= 0] <- 0
  out[C >= 1] <- 1
  out[!is.finite(C)] <- NA_real_
  dimnames(out) <- dimnames(C)
  attr(out, "integration_formula") <- paste(
    "C_multiome = C_RNA * 2^R / (1 - C_RNA + C_RNA * 2^R);",
    "non-finite Pando R := 0"
  )
  attr(out, "regulatory_alpha") <- 1
  attr(out, "rna_only_fallback_mask") <- fallback
  attr(out, "fallback_policy") <-
    "nonfinite_pando_modifier_uses_neutral_R_and_equals_RNA_only"
  out
}

.rc_cell_first_projection_layer1_base <- .rc_cell_first_projection_layer1
.rc_cell_first_projection_layer1 <- function(
    grn_result, metacell_object, membership, metacell_meta, gem,
    condition_col, celltype_col, rna_assay,
    projection_component, comparison_support, regulatory_alpha = 1,
    gpr_and_method, gene_half_saturation, parallel, BPPARAM) {
  if (!is.numeric(regulatory_alpha) || length(regulatory_alpha) != 1L ||
      !is.finite(regulatory_alpha) ||
      !isTRUE(all.equal(as.numeric(regulatory_alpha), 1))) {
    stop("Canonical RegCompass requires `regulatory_alpha = 1`.",
         call. = FALSE)
  }
  id_col <- if ("metacell_id" %in% colnames(metacell_meta)) {
    "metacell_id"
  } else if ("pool_id" %in% colnames(metacell_meta)) {
    "pool_id"
  } else {
    stop("Metacell metadata lack metacell_id/pool_id.", call. = FALSE)
  }
  unit_id <- as.character(metacell_meta[[id_col]])
  cell_type <- stats::setNames(
    as.character(metacell_meta[[celltype_col]]), unit_id
  )
  old_cell_type <- .rc_layer1_policy_state$cell_type
  .rc_layer1_policy_state$cell_type <- cell_type
  on.exit({
    if (is.null(old_cell_type)) {
      rm(list = "cell_type", envir = .rc_layer1_policy_state)
    } else {
      .rc_layer1_policy_state$cell_type <- old_cell_type
    }
  }, add = TRUE)

  out <- .rc_cell_first_projection_layer1_base(
    grn_result = grn_result,
    metacell_object = metacell_object,
    membership = membership,
    metacell_meta = metacell_meta,
    gem = gem,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    projection_component = projection_component,
    comparison_support = comparison_support,
    regulatory_alpha = 1,
    gpr_and_method = gpr_and_method,
    gene_half_saturation = gene_half_saturation,
    parallel = parallel,
    BPPARAM = BPPARAM
  )
  modifier <- as.matrix(out$gene_regulatory_modifier_common_oof)
  fallback <- !is.finite(modifier)
  out$regulatory_fallback <- list(
    policy = "rna_only_for_nonfinite_pando_modifier",
    neutral_modifier = 0,
    regulatory_alpha = 1,
    gene_metacell_mask = fallback,
    n_fallback = sum(fallback),
    fallback_fraction = mean(fallback),
    result_annotation = paste(
      "Pando modifier unavailable: support and reaction expression use",
      "the RNA-only value for that gene-metacell entry"
    )
  )
  out$capacity_params$regulatory_alpha <- 1
  out$capacity_params$regulatory_odds_budget <- 2
  out$capacity_params$or_method <- "sum"
  out$capacity_params$or_method_semantics <-
    "COMPASS isoform summing across complete OR branches"
  out$projection_provenance$unavailable_target_policy <-
    "rna_only_neutral_modifier_fallback"
  out$projection_provenance$rna_only_fallback_reported <- TRUE
  out$projection_provenance$reference_condition_contrast <- FALSE
  out$depth_diagnostics$latent_expression_model <-
    "gamma_poisson_empirical_bayes_by_cell_type_v2"
  out$depth_diagnostics$prior_estimation_scope <- "gene_by_cell_type"
  out
}

.rc_regcompass_step_layer1_base <- rc_regcompass_step_layer1
rc_regcompass_step_layer1 <- function(
    grn, metacells, meta_modules, gem, outdir,
    projection_component = "condition",
    comparison_support = c("auto", "pairwise_common", "global_common"),
    regulatory_alpha = 1,
    gpr_and_method = c("min", "median", "mean"),
    gene_half_saturation = getOption("RegCompassR.cpm_half_saturation", 1),
    parallel = TRUE,
    BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  if (!is.numeric(regulatory_alpha) || length(regulatory_alpha) != 1L ||
      !is.finite(regulatory_alpha) ||
      !isTRUE(all.equal(as.numeric(regulatory_alpha), 1))) {
    stop("Canonical RegCompass requires `regulatory_alpha = 1`.",
         call. = FALSE)
  }
  .rc_regcompass_step_layer1_base(
    grn = grn,
    metacells = metacells,
    meta_modules = meta_modules,
    gem = gem,
    outdir = outdir,
    projection_component = projection_component,
    comparison_support = comparison_support,
    regulatory_alpha = 1,
    gpr_and_method = gpr_and_method,
    gene_half_saturation = gene_half_saturation,
    parallel = parallel,
    BPPARAM = BPPARAM,
    progress = progress
  )
}
