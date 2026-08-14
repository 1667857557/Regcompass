# Layer 1 RNA and regulatory support integration.

.rc_latent_metacell_expression <- function(
    counts, library_size, mu_min = 0.1, cell_type) {
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
    stop("Latent RNA estimation requires raw integer-like counts.",
         call. = FALSE)
  }
  if (!is.null(names(cell_type))) {
    cell_type <- as.character(cell_type[colnames(counts)])
  } else {
    cell_type <- as.character(cell_type)
  }
  if (length(cell_type) != ncol(counts) || anyNA(cell_type)) {
    stop("Cell-type labels do not align to metacell RNA counts.",
         call. = FALSE)
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
  posterior_shape <- observed_cpm + shape_unit
  posterior_rate <- rate_unit + 1
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
  observed_zero <- counts == 0
  zero_class <- matrix(
    "observed_positive", nrow(counts), ncol(counts),
    dimnames = dimnames(counts)
  )
  zero_class[observed_zero] <- "observed_zero_continuous_eb"
  observation_weight <- 1 / posterior_rate
  prior_weight <- rate_unit / posterior_rate
  list(
    latent_cpm = latent_cpm,
    latent_log_expression = log1p(latent_cpm),
    posterior_positive_probability = positive,
    posterior_zero_probability = 1 - positive,
    zero_class = zero_class,
    observed_zero = observed_zero,
    prior_mean = prior_mean,
    prior_variance = prior_variance,
    prior_shape = prior_shape,
    prior_rate = prior_rate,
    prior_weight = prior_weight,
    observation_weight = observation_weight,
    model = "normalized_unit_gamma_empirical_bayes_by_cell_type_v4",
    prior_estimation_scope = "gene_by_cell_type",
    posterior_update_scope =
      "one_normalized_metacell_unit_independent_of_library_size"
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
    stop("Canonical RegCompass requires `regulatory_alpha = 1`.",
         call. = FALSE)
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

.rc_integrate_regulatory_expression <- function(
    rna_expression, regulatory_modifier) {
  rna_expression <- as.matrix(rna_expression)
  regulatory_modifier <- as.matrix(regulatory_modifier)
  if (!identical(dimnames(rna_expression), dimnames(regulatory_modifier))) {
    stop(
      "Quantitative RNA expression and regulatory modifier matrices must align exactly.",
      call. = FALSE
    )
  }
  observed <- is.finite(rna_expression)
  X <- pmax(rna_expression, 0)
  fallback <- !is.finite(regulatory_modifier)
  R <- regulatory_modifier
  R[fallback] <- 0
  R <- pmin(pmax(R, -1), 1)
  out <- X * 2^R
  out[observed & X <= 0] <- 0
  out[!observed] <- NA_real_
  dimnames(out) <- dimnames(X)
  attr(out, "rna_only_fallback_mask") <- fallback
  attr(out, "integration_formula") <-
    "X_multiome=X_RNA*2^R; nonfinite R:=0; R clipped to [-1,1]"
  out
}

.rc_projection_scale <- function(projection, unit_meta, celltype_col) {
  scale <- matrix(
    NA_real_, nrow(projection), ncol(projection),
    dimnames = dimnames(projection)
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
          sqrt(mean(value^2)),
          1e-6,
          na.rm = TRUE
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

.rc_scaled_regulatory_modifier <- function(projection, reliability, scale) {
  if (!identical(dimnames(projection), dimnames(reliability)) ||
      !identical(dimnames(projection), dimnames(scale))) {
    stop("Projection, reliability, and calibration scale must align.",
         call. = FALSE)
  }
  projection <- as.matrix(projection)
  reliability <- as.matrix(reliability)
  scale <- as.matrix(scale)
  value <- matrix(
    NA_real_, nrow(projection), ncol(projection), dimnames = dimnames(projection)
  )
  available <- is.finite(reliability)
  neutral <- available & reliability == 0
  value[neutral] <- 0
  regulated <- available & reliability != 0 &
    is.finite(projection) & is.finite(scale) & scale > 0
  value[regulated] <- reliability[regulated] *
    tanh(projection[regulated] / scale[regulated])
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
