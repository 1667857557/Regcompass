.rc_latent_metacell_expression <- function(
    counts, library_size, mu_min = 0.1,
    structural_zero_probability = 0.05) {
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
  integer_tolerance <- sqrt(.Machine$double.eps) *
    pmax(1, abs(counts))
  if (any(abs(counts - round(counts)) > integer_tolerance)) {
    stop(
      "Gamma-Poisson latent RNA estimation requires raw integer-like counts, ",
      "not normalized expression.",
      call. = FALSE
    )
  }
  exposure <- library_size / 1e6
  observed_cpm <- sweep(counts, 2L, exposure, "/")
  prior_mean <- rowMeans(observed_cpm)
  prior_variance <- apply(observed_cpm, 1L, stats::var)
  prior_variance[!is.finite(prior_variance) | prior_variance <= 0] <-
    pmax(prior_mean[
      !is.finite(prior_variance) | prior_variance <= 0
    ]^2 / 1e6, 1e-12)
  prior_shape <- pmin(
    1e6, pmax(0.1, prior_mean^2 / pmax(prior_variance, 1e-12))
  )
  prior_rate <- prior_shape / pmax(prior_mean, 1e-8)
  posterior_shape <- sweep(counts, 1L, prior_shape, "+")
  posterior_rate <- outer(prior_rate, exposure, "+")
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
  zero_class[observed_zero & !credible_structural] <- "sampling_limited_zero"
  zero_class[credible_structural] <- "credible_structural_zero"
  latent_cpm[credible_structural] <- 0
  dimnames(latent_cpm) <- dimnames(counts)
  latent_log <- log1p(latent_cpm)
  list(
    latent_cpm = latent_cpm,
    latent_log_expression = latent_log,
    posterior_positive_probability = posterior_positive,
    zero_class = zero_class,
    observed_zero = observed_zero,
    dispersion = stats::setNames(1 / prior_shape, rownames(counts)),
    prior_shape = stats::setNames(prior_shape, rownames(counts)),
    prior_rate = stats::setNames(prior_rate, rownames(counts)),
    model = "gamma_poisson_empirical_bayes_with_library_offset_v1",
    marginal_count_model = "negative_binomial",
    mu_min = mu_min,
    structural_zero_probability = structural_zero_probability
  )
}

.rc_depth_matched_counts <- function(
    counts, library_size, unit_meta, celltype_col, seed = 12345L) {
  counts <- as.matrix(counts)
  celltype <- as.character(unit_meta[[celltype_col]])
  target_depth <- vapply(unique(celltype), function(value) {
    min(library_size[celltype == value])
  }, numeric(1))
  names(target_depth) <- unique(celltype)
  probability <- pmin(1, target_depth[celltype] / library_size)
  set.seed(seed)
  sampled <- matrix(
    stats::rbinom(
      length(counts),
      size = as.integer(round(as.numeric(counts))),
      prob = rep(probability, each = nrow(counts))
    ),
    nrow = nrow(counts),
    dimnames = dimnames(counts)
  )
  list(
    counts = sampled,
    library_size = stats::setNames(
      target_depth[celltype], colnames(counts)
    ),
    target_depth = target_depth
  )
}

.rc_common_depth_interval <- function(
    library_size, unit_meta, condition_col, celltype_col) {
  condition <- as.character(unit_meta[[condition_col]])
  celltype <- as.character(unit_meta[[celltype_col]])
  keep <- rep(FALSE, length(library_size))
  diagnostics <- list()
  for (value in unique(celltype)) {
    ct <- celltype == value
    ranges <- lapply(unique(condition[ct]), function(level) {
      depth <- library_size[ct & condition == level]
      stats::quantile(depth, c(0.1, 0.9), names = FALSE, na.rm = TRUE)
    })
    lower <- max(vapply(ranges, `[[`, numeric(1), 1L))
    upper <- min(vapply(ranges, `[[`, numeric(1), 2L))
    overlap <- is.finite(lower) && is.finite(upper) && lower <= upper
    if (overlap) {
      keep[ct] <- library_size[ct] >= lower & library_size[ct] <= upper
    }
    diagnostics[[length(diagnostics) + 1L]] <- data.frame(
      cell_type = value,
      lower_library_size = lower,
      upper_library_size = upper,
      condition_interval_overlap = overlap,
      fallback_used = FALSE,
      retained_metacells = sum(keep[ct]),
      total_metacells = sum(ct),
      stringsAsFactors = FALSE
    )
  }
  list(keep = keep, diagnostics = do.call(rbind, diagnostics))
}

.rc_depth_direction_diagnostics <- function(
    full, downsampled, interval, unit_meta, condition_col, celltype_col) {
  condition <- as.character(unit_meta[[condition_col]])
  celltype <- as.character(unit_meta[[celltype_col]])
  rows <- list()
  for (value in unique(celltype)) {
    selected_condition <- unique(condition[celltype == value])
    if (length(selected_condition) != 2L) next
    a <- celltype == value & condition == selected_condition[[1L]]
    b <- celltype == value & condition == selected_condition[[2L]]
    contrast <- function(x) {
      rowMeans(x[, b, drop = FALSE], na.rm = TRUE) -
        rowMeans(x[, a, drop = FALSE], na.rm = TRUE)
    }
    full_effect <- contrast(full)
    downsampled_effect <- contrast(downsampled)
    interval_effect <- contrast(interval)
    signs <- cbind(
      full = sign(full_effect),
      downsampled = sign(downsampled_effect),
      common_interval = sign(interval_effect)
    )
    finite_count <- rowSums(is.finite(signs))
    stable <- apply(signs, 1L, function(x) {
      x <- x[is.finite(x)]
      length(x) >= 2L && length(unique(x)) == 1L
    })
    rows[[length(rows) + 1L]] <- data.frame(
      reaction_id = rownames(full),
      cell_type = value,
      condition_a = selected_condition[[1L]],
      condition_b = selected_condition[[2L]],
      full_depth_effect = full_effect,
      downsampled_effect = downsampled_effect,
      common_interval_effect = interval_effect,
      finite_sensitivity_routes = finite_count,
      depth_sensitivity_flag = finite_count >= 2L & !stable,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

.rc_projection_scale <- function(
    projection_common, unit_meta, celltype_col) {
  celltypes <- unique(as.character(unit_meta[[celltype_col]]))
  scale <- matrix(
    NA_real_, nrow(projection_common), ncol(projection_common),
    dimnames = dimnames(projection_common)
  )
  diagnostics <- list()
  for (gene in rownames(projection_common)) {
    for (celltype in celltypes) {
      unit <- as.character(unit_meta$unit_id[
        as.character(unit_meta[[celltype_col]]) == celltype
      ])
      primary <- projection_common[gene, unit]
      value <- primary[is.finite(primary)]
      source <- "common_oof"
      if (length(value) < 2L || diff(range(value)) <= 1e-12) {
        value <- numeric()
        source <- "fixed_unit_common_oof_fallback"
      }
      robust_scale <- if (length(value) >= 2L) {
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
      scale[gene, unit] <- robust_scale
      diagnostics[[length(diagnostics) + 1L]] <- data.frame(
        target = gene,
        cell_type = celltype,
        projection_scale = robust_scale,
        scale_source = source,
        projection_q01 = if (length(value)) {
          unname(stats::quantile(value, 0.01))
        } else NA_real_,
        projection_q50 = if (length(value)) stats::median(value) else NA_real_,
        projection_q99 = if (length(value)) {
          unname(stats::quantile(value, 0.99))
        } else NA_real_,
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
  modifier <- reliability * tanh(projection / scale)
  modifier[!is.finite(projection) | !is.finite(reliability)] <- NA_real_
  dimnames(modifier) <- dimnames(projection)
  modifier
}

.rc_reaction_support_contract <- function(
    parsed_gpr, gene_support, reaction_expression, pathway) {
  genes_available <- rownames(gene_support)
  reaction_ids <- names(parsed_gpr)
  units <- colnames(gene_support)
  rows <- lapply(reaction_ids, function(reaction_id) {
    groups <- parsed_gpr[[reaction_id]]
    required <- unique(tolower(unlist(groups, use.names = FALSE)))
    complete_group <- vapply(groups, function(group) {
      all(unique(tolower(group)) %in% genes_available)
    }, logical(1))
    has_complete_isozyme <- any(complete_group)
    partial_isozyme_coverage <- has_complete_isozyme &&
      any(!complete_group)
    do.call(rbind, lapply(units, function(unit) {
      value <- reaction_expression[reaction_id, unit]
      reason <- if (!length(groups) || !length(required)) {
        "no_gpr"
      } else if (!has_complete_isozyme) {
        "gpr_mapping_incomplete"
      } else if (!is.finite(value)) {
        if (identical(pathway, "rna_only")) {
          "rna_gene_missing"
        } else {
          "grn_unavailable"
        }
      } else if (value <= 0) {
        "observed_zero"
      } else {
        "observed_positive"
      }
      data.frame(
        reaction_id = reaction_id,
        metacell_id = unit,
        pathway = pathway,
        reaction_expression = value,
        reaction_expression_available = is.finite(value),
        complete_isozyme_groups = sum(complete_group),
        total_isozyme_groups = length(groups),
        partial_isozyme_coverage = partial_isozyme_coverage,
        reaction_missing_reason = if (reason %in% c(
          "observed_zero", "observed_positive"
        )) NA_character_ else reason,
        reaction_support_class = reason,
        stringsAsFactors = FALSE
      )
    }))
  })
  do.call(rbind, rows)
}

.rc_gpr_best_group_fraction <- function(parsed_gpr, gene_available) {
  gene_available <- as.matrix(gene_available)
  if (!is.logical(gene_available) || is.null(rownames(gene_available)) ||
      is.null(colnames(gene_available))) {
    stop("GPR availability requires a named logical gene-by-unit matrix.",
         call. = FALSE)
  }
  gene_names <- tolower(rownames(gene_available))
  answer <- matrix(
    NA_real_, length(parsed_gpr), ncol(gene_available),
    dimnames = list(names(parsed_gpr), colnames(gene_available))
  )
  for (reaction_id in names(parsed_gpr)) {
    groups <- parsed_gpr[[reaction_id]]
    if (!length(groups)) next
    group_fraction <- vapply(groups, function(group) {
      requested <- unique(tolower(group))
      index <- match(requested, gene_names)
      vapply(seq_len(ncol(gene_available)), function(unit) {
        observed <- rep(FALSE, length(requested))
        mapped <- !is.na(index)
        observed[mapped] <- gene_available[index[mapped], unit]
        mean(observed)
      }, numeric(1))
    }, numeric(ncol(gene_available)))
    if (is.null(dim(group_fraction))) {
      group_fraction <- matrix(group_fraction, nrow = ncol(gene_available))
    }
    answer[reaction_id, ] <- apply(group_fraction, 1L, max)
  }
  answer
}

.rc_zero_pattern_diagnostics <- function(
    counts, zero_class, unit_meta, condition_col, celltype_col,
    asymmetry_threshold = 0.25) {
  condition <- as.character(unit_meta[[condition_col]])
  celltype <- as.character(unit_meta[[celltype_col]])
  rows <- list()
  for (gene in rownames(counts)) {
    for (value in unique(celltype)) {
      selected <- celltype == value
      levels_condition <- unique(condition[selected])
      observed_fraction <- vapply(levels_condition, function(level) {
        index <- selected & condition == level
        mean(counts[gene, index] == 0)
      }, numeric(1))
      sampling_fraction <- vapply(levels_condition, function(level) {
        index <- selected & condition == level
        mean(zero_class[gene, index] == "sampling_limited_zero")
      }, numeric(1))
      structural_fraction <- vapply(levels_condition, function(level) {
        index <- selected & condition == level
        mean(zero_class[gene, index] == "credible_structural_zero")
      }, numeric(1))
      asymmetry <- if (length(observed_fraction) >= 2L) {
        diff(range(observed_fraction))
      } else {
        NA_real_
      }
      rows[[length(rows) + 1L]] <- data.frame(
        target = gene,
        cell_type = value,
        observed_zero_fraction = mean(counts[gene, selected] == 0),
        posterior_sampling_zero_fraction = mean(
          zero_class[gene, selected] == "sampling_limited_zero"
        ),
        credible_structural_zero_fraction = mean(
          zero_class[gene, selected] == "credible_structural_zero"
        ),
        always_zero_condition_fraction = mean(observed_fraction == 1),
        zero_pattern_asymmetry = asymmetry,
        missing_target_fraction = mean(!is.finite(counts[gene, selected])),
        zero_support_sensitive =
          is.finite(asymmetry) && asymmetry > asymmetry_threshold,
        zero_asymmetry_threshold = asymmetry_threshold,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

.rc_reaction_sensitivity_from_genes <- function(
    parsed_gpr, gene_diagnostics, flag_col) {
  if (!is.data.frame(gene_diagnostics) ||
      !all(c("target", "cell_type", flag_col) %in%
           colnames(gene_diagnostics))) {
    stop("Gene sensitivity diagnostics are incomplete.", call. = FALSE)
  }
  rows <- list()
  for (reaction_id in names(parsed_gpr)) {
    genes <- unique(tolower(unlist(
      parsed_gpr[[reaction_id]], use.names = FALSE
    )))
    for (celltype in unique(as.character(gene_diagnostics$cell_type))) {
      selected <- tolower(as.character(gene_diagnostics$target)) %in% genes &
        as.character(gene_diagnostics$cell_type) == celltype
      rows[[length(rows) + 1L]] <- data.frame(
        reaction_id = reaction_id,
        cell_type = celltype,
        sensitivity_flag = any(
          gene_diagnostics[[flag_col]][selected] %in% TRUE
        ),
        contributing_target_count = sum(selected),
        stringsAsFactors = FALSE
      )
    }
  }
  answer <- do.call(rbind, rows)
  colnames(answer)[colnames(answer) == "sensitivity_flag"] <- flag_col
  answer
}

.rc_cell_first_projection_layer1 <- function(
    grn_result, metacell_object, membership, metacell_meta, gem,
    condition_col, celltype_col, rna_assay,
    projection_component, comparison_support, regulatory_alpha,
    gpr_and_method, gene_half_saturation, parallel, BPPARAM) {
  if (!identical(
    grn_result$schema_version, "regcompass_condition_grn_fit_v5"
  ) || !inherits(grn_result$pando_grn_data, "GRNData")) {
    stop("Stage 1 lacks the exact Pando ConditionGRNFit v5 artifact.",
         call. = FALSE)
  }
  if (!identical(projection_component, "condition")) {
    stop("The canonical penalty path requires projection_component = 'condition'.",
         call. = FALSE)
  }
  if (comparison_support %in% c("condition_estimable", "strict")) {
    stop("Diagnostic support policies cannot enter the primary penalty path.",
         call. = FALSE)
  }
  if (!is.data.frame(membership) ||
      !all(c("cell_id", "metacell_id") %in% colnames(membership)) ||
      anyDuplicated(as.character(membership$cell_id))) {
    stop("SuperCell membership must map every cell exactly once.",
         call. = FALSE)
  }
  units <- as.character(colnames(metacell_object))
  if (!setequal(unique(as.character(membership$metacell_id)), units)) {
    stop("SuperCell membership and metacell assay IDs differ.", call. = FALSE)
  }
  id_col <- if ("metacell_id" %in% colnames(metacell_meta)) {
    "metacell_id"
  } else if ("pool_id" %in% colnames(metacell_meta)) {
    "pool_id"
  } else {
    stop("SuperCell metadata lack metacell_id/pool_id.", call. = FALSE)
  }
  unit_meta <- as.data.frame(metacell_meta)
  unit_meta$unit_id <- as.character(unit_meta[[id_col]])
  if (!all(c(condition_col, celltype_col) %in% colnames(unit_meta)) ||
      anyDuplicated(unit_meta$unit_id) ||
      !setequal(unit_meta$unit_id, units)) {
    stop("SuperCell metadata do not define aligned condition/cell-type units.",
         call. = FALSE)
  }
  unit_meta <- unit_meta[match(units, unit_meta$unit_id), , drop = FALSE]
  unit_meta$pool_id <- unit_meta$unit_id
  membership_count <- table(as.character(membership$metacell_id))
  unit_meta$metacell_cell_count <- as.integer(membership_count[units])

  parsed <- rc_parse_gpr_table(gem$gpr_table)
  gpr_genes <- unique(tolower(unlist(parsed, use.names = FALSE)))
  counts <- .rc_get_assay_counts(metacell_object, rna_assay)
  full_library_size <- Matrix::colSums(counts)
  rna_counts <- counts[
    tolower(rownames(counts)) %in% gpr_genes, units, drop = FALSE
  ]
  rownames(rna_counts) <- tolower(rownames(rna_counts))
  if (anyDuplicated(rownames(rna_counts))) {
    stop("Duplicated GPR genes after case normalization.", call. = FALSE)
  }
  latent <- .rc_latent_metacell_expression(
    rna_counts, full_library_size[units]
  )
  depth_matched <- .rc_depth_matched_counts(
    rna_counts, full_library_size[units], unit_meta, celltype_col
  )
  latent_depth_matched <- .rc_latent_metacell_expression(
    depth_matched$counts, depth_matched$library_size
  )
  genes <- rownames(latent$latent_log_expression)
  empty <- matrix(
    NA_real_, length(genes), length(units),
    dimnames = list(genes, units)
  )
  projection_common <- projection_full <- reliability <- empty
  reliability_available <- matrix(
    FALSE, length(genes), length(units), dimnames = dimnames(empty)
  )
  coverage <- list()
  fits <- grn_result$condition_grn_fits
  if (!is.list(fits) || !length(fits)) {
    stop("Stage 1 contains no Pando ConditionGRNFit v5 contracts.",
         call. = FALSE)
  }
  assign_projection <- function(
      destination, projected, primary_penalty = TRUE) {
    if (!identical(
      projected$projection_origin,
      "outer_condition_stratified_cell_oof"
    ) || isTRUE(projected$full_fit_projection_used_for_penalty)) {
      stop("Regcompass accepts only outer-heldout Pando projections.",
           call. = FALSE)
    }
    if (!identical(
          isTRUE(projected$projection_used_for_penalty),
          isTRUE(primary_penalty)
        )) {
      stop(
        if (isTRUE(primary_penalty)) {
          "The common-support projection is not eligible for the primary penalty."
        } else {
          "The condition-full comparator was incorrectly marked as a primary penalty."
        },
        call. = FALSE
      )
    }
    score <- t(as.matrix(projected$gene_score))
    rownames(score) <- tolower(rownames(score))
    target <- intersect(rownames(score), rownames(destination))
    group <- intersect(colnames(score), colnames(destination))
    if (!setequal(colnames(score), group)) {
      stop("Pando group projections contain unknown metacells.",
           call. = FALSE)
    }
    destination[target, group] <- score[target, group, drop = FALSE]
    destination
  }
  for (fit in fits) {
    if (!identical(fit$schema_version, "pando_condition_grn_fit_v5") ||
        !identical(
          fit$projection_origin, "outer_condition_stratified_cell_oof"
        )) {
      stop("Layer 1 requires ConditionGRNFit v5 OOF contracts.",
           call. = FALSE)
    }
    resolved <- if (identical(comparison_support, "auto")) {
      if (length(fit$condition_levels) == 2L) {
        "pairwise_common"
      } else {
        "global_common"
      }
    } else {
      comparison_support
    }
    comparison_conditions <- if (resolved == "pairwise_common") {
      fit$comparison_conditions
    } else {
      NULL
    }
    project <- function(support, diagnostic = FALSE) {
      cell_projection <- Pando::project_condition_grn_cells(
        object = grn_result$pando_grn_data,
        fit = fit,
        component = "condition",
        scale = "std",
        support_policy = support,
        comparison_conditions = comparison_conditions,
        origin = "oof",
        diagnostic_only = diagnostic
      )
      Pando::aggregate_condition_grn_projection(
        cell_projection, membership, group_col = "metacell_id"
      )
    }
    common <- project(resolved)
    full <- project("condition_estimable", diagnostic = TRUE)
    projection_common <- assign_projection(
      projection_common, common, primary_penalty = TRUE
    )
    projection_full <- assign_projection(
      projection_full, full, primary_penalty = FALSE
    )
    fit_units <- rownames(common$gene_score)
    expected <- as.character(unit_meta$unit_id[
      as.character(unit_meta[[celltype_col]]) == fit$cell_type &
        as.character(unit_meta[[condition_col]]) %in% fit$condition_levels
    ])
    if (!setequal(fit_units, expected)) {
      stop("Pando paired cells do not cover every fitted SuperCell.",
           call. = FALSE)
    }
    pooled_oof <- fit$target_rsq_oof_pooled
    names(pooled_oof) <- tolower(names(pooled_oof))
    target <- intersect(names(pooled_oof), genes)
    available <- is.finite(pooled_oof[target]) &
      fit$predictive_oof_available[
        match(target, tolower(names(fit$predictive_oof_available)))
      ] &
      fit$oof_cell_coverage[
        match(target, tolower(names(fit$oof_cell_coverage)))
      ] == 1
    available[is.na(available)] <- FALSE
    q <- sqrt(pmin(1, pmax(0, as.numeric(pooled_oof[target]))))
    q[!available] <- NA_real_
    reliability[target, fit_units] <- q
    reliability_available[target, fit_units] <- available
    status <- common$source_projection$target_condition_status
    status$comparison_support <- resolved
    status$projection_origin <- common$projection_origin
    status$projection_used_for_penalty <-
      common$projection_used_for_penalty
    status$oof_cell_coverage <- fit$oof_cell_coverage[
      match(
        tolower(status$target),
        tolower(names(fit$oof_cell_coverage))
      )
    ]
    status$oof_projection_available_fraction <-
      fit$oof_projection_available_fraction[
        match(
          tolower(status$target),
          tolower(names(fit$oof_projection_available_fraction))
        )
      ]
    coverage[[length(coverage) + 1L]] <- status
  }

  calibration <- .rc_projection_scale(
    projection_common, unit_meta, celltype_col
  )
  modifier_common <- .rc_scaled_oof_modifier(
    projection_common, reliability, calibration$scale
  )
  modifier_full <- .rc_scaled_oof_modifier(
    projection_full, reliability, calibration$scale
  )
  calibration$diagnostics$modifier_q01 <- NA_real_
  calibration$diagnostics$modifier_q50 <- NA_real_
  calibration$diagnostics$modifier_q99 <- NA_real_
  for (i in seq_len(nrow(calibration$diagnostics))) {
    gene <- calibration$diagnostics$target[[i]]
    celltype <- calibration$diagnostics$cell_type[[i]]
    selected <- as.character(unit_meta$unit_id[
      as.character(unit_meta[[celltype_col]]) == celltype
    ])
    value <- modifier_common[gene, selected]
    value <- value[is.finite(value)]
    if (length(value)) {
      quantile <- stats::quantile(
        value, c(0.01, 0.5, 0.99), names = FALSE
      )
      calibration$diagnostics$modifier_q01[[i]] <- quantile[[1L]]
      calibration$diagnostics$modifier_q50[[i]] <- quantile[[2L]]
      calibration$diagnostics$modifier_q99[[i]] <- quantile[[3L]]
    }
  }
  gene_support_rna <- rc_gene_score(
    latent$latent_log_expression,
    mode = "absolute",
    half_saturation = gene_half_saturation
  )
  gene_support_rna_depth_matched <- rc_gene_score(
    latent_depth_matched$latent_log_expression,
    mode = "absolute",
    half_saturation = gene_half_saturation
  )
  gene_support_common <- .rc_integrate_regulatory_support(
    gene_support_rna, modifier_common, alpha = regulatory_alpha
  )
  gene_support_full <- .rc_integrate_regulatory_support(
    gene_support_rna, modifier_full, alpha = regulatory_alpha
  )
  reaction_rna <- rc_reaction_capacity(
    parsed, gene_support_rna, promiscuity_mode = "none",
    and_method = gpr_and_method, or_method = "sum",
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )
  reaction_common <- rc_reaction_capacity(
    parsed, gene_support_common, promiscuity_mode = "none",
    and_method = gpr_and_method, or_method = "sum",
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )
  reaction_full <- rc_reaction_capacity(
    parsed, gene_support_full, promiscuity_mode = "none",
    and_method = gpr_and_method, or_method = "sum",
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )
  reaction_depth_matched <- rc_reaction_capacity(
    parsed, gene_support_rna_depth_matched, promiscuity_mode = "none",
    and_method = gpr_and_method, or_method = "sum",
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )
  alpha_grid <- unique(c(0, 0.25, 0.5, 1, regulatory_alpha))
  alpha_sensitivity <- stats::setNames(
    lapply(alpha_grid, function(alpha_value) {
      support <- .rc_integrate_regulatory_support(
        gene_support_rna, modifier_common, alpha = alpha_value
      )
      rc_reaction_capacity(
        parsed, support, promiscuity_mode = "none",
        and_method = gpr_and_method, or_method = "sum",
        BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
      )
    }),
    format(alpha_grid, trim = TRUE, scientific = FALSE)
  )
  depth_interval <- .rc_common_depth_interval(
    full_library_size[units], unit_meta, condition_col, celltype_col
  )
  reaction_common_depth_interval <- reaction_rna
  reaction_common_depth_interval[, !depth_interval$keep] <- NA_real_
  depth_direction <- .rc_depth_direction_diagnostics(
    reaction_rna,
    reaction_depth_matched,
    reaction_common_depth_interval,
    unit_meta,
    condition_col,
    celltype_col
  )
  support_contract <- do.call(rbind, list(
    .rc_reaction_support_contract(
      parsed, gene_support_rna, reaction_rna, "rna_only"
    ),
    .rc_reaction_support_contract(
      parsed, gene_support_common, reaction_common, "common_oof"
    ),
    .rc_reaction_support_contract(
      parsed, gene_support_full, reaction_full, "condition_full_oof"
    )
  ))
  unit_meta$rna_library_size <- as.numeric(full_library_size[units])
  fragment_col <- intersect(
    c("atac_fragment_size", "atac_fragments", "total_fragments"),
    colnames(unit_meta)
  )
  unit_meta$atac_fragment_size <- if (length(fragment_col)) {
    as.numeric(unit_meta[[fragment_col[[1L]]]])
  } else {
    NA_real_
  }
  calibration$diagnostics$link_function <- "tanh(G / shared_oof_scale)"
  calibration$diagnostics$link_saturation_sensitive <- vapply(
    seq_len(nrow(calibration$diagnostics)), function(i) {
      gene <- calibration$diagnostics$target[[i]]
      celltype <- calibration$diagnostics$cell_type[[i]]
      selected <- as.character(unit_meta$unit_id[
        as.character(unit_meta[[celltype_col]]) == celltype
      ])
      value <- modifier_common[gene, selected]
      value <- value[is.finite(value)]
      length(value) > 0L && mean(abs(value) > 0.95) > 0.01
    }, logical(1)
  )
  zero_pattern <- .rc_zero_pattern_diagnostics(
    rna_counts, latent$zero_class, unit_meta, condition_col, celltype_col
  )
  common_support_fraction <- .rc_gpr_best_group_fraction(
    parsed, is.finite(modifier_common)
  )
  condition_full_support_fraction <- .rc_gpr_best_group_fraction(
    parsed, is.finite(modifier_full)
  )
  reaction_zero_sensitivity <- .rc_reaction_sensitivity_from_genes(
    parsed, zero_pattern, "zero_support_sensitive"
  )
  reaction_link_sensitivity <- .rc_reaction_sensitivity_from_genes(
    parsed, calibration$diagnostics, "link_saturation_sensitive"
  )
  list(
    schema_version = "regcompass_condition_grn_layer1_v4",
    reaction_expression = reaction_common,
    reaction_expression_common_oof = reaction_common,
    reaction_expression_condition_full_oof = reaction_full,
    reaction_expression_rna_only = reaction_rna,
    reaction_expression_depth_matched_rna = reaction_depth_matched,
    reaction_expression_common_depth_interval_rna =
      reaction_common_depth_interval,
    reaction_expression_alpha_sensitivity = alpha_sensitivity,
    reaction_support_contract = support_contract,
    reaction_expression_available = is.finite(reaction_common),
    rna_metacell_latent_log_expression =
      latent$latent_log_expression,
    rna_metacell_latent_cpm = latent$latent_cpm,
    posterior_positive_probability =
      latent$posterior_positive_probability,
    posterior_zero_probability =
      1 - latent$posterior_positive_probability,
    rna_zero_class = latent$zero_class,
    gene_support_rna = gene_support_rna,
    gene_support_common_oof = gene_support_common,
    gene_support_condition_full_oof = gene_support_full,
    gene_projection_common_oof = projection_common,
    gene_projection_condition_full_oof = projection_full,
    gene_projection_raw = projection_common,
    gene_projection_scale = calibration$scale,
    gene_regulatory_reliability = reliability,
    gene_regulatory_reliability_available = reliability_available,
    gene_regulatory_available = is.finite(projection_common),
    gene_regulatory_modifier = modifier_common,
    gene_regulatory_modifier_common_oof = modifier_common,
    gene_regulatory_modifier_condition_full_oof = modifier_full,
    gene_support_multiome = gene_support_common,
    projection_coverage = do.call(rbind, coverage),
    projection_calibration = calibration$diagnostics,
    reaction_common_support_fraction = common_support_fraction,
    reaction_condition_full_support_fraction =
      condition_full_support_fraction,
    zero_pattern_diagnostics = zero_pattern,
    reaction_zero_support_sensitivity = reaction_zero_sensitivity,
    reaction_link_saturation_sensitivity = reaction_link_sensitivity,
    parsed_gpr = parsed,
    gpr_diagnostics = rc_gpr_diagnostics(parsed, genes),
    unit_meta = unit_meta,
    metacell_meta = unit_meta,
    layer1_unit = "SuperCell_condition_by_broad_cell_type",
    depth_diagnostics = list(
      rna_library_size = stats::setNames(
        as.numeric(full_library_size[units]), units
      ),
      atac_fragment_size = stats::setNames(
        unit_meta$atac_fragment_size, units
      ),
      metacell_cell_count = stats::setNames(
        unit_meta$metacell_cell_count, units
      ),
      depth_balance_policy =
        "shared_celltype_targets_plus_full_depth_and_downsample_sensitivity",
      latent_expression_model = latent$model,
      dispersion_model = "gene_specific_gamma_poisson_empirical_bayes",
      depth_matched_library_size = depth_matched$library_size,
      common_depth_interval = depth_interval$diagnostics,
      common_depth_interval_keep = stats::setNames(
        depth_interval$keep, units
      ),
      reaction_depth_sensitivity = depth_direction
    ),
    zero_diagnostics = list(
      observed_zero_fraction = rowMeans(latent$observed_zero),
      posterior_sampling_zero_fraction = rowMeans(
        latent$zero_class == "sampling_limited_zero"
      ),
      credible_structural_zero_fraction = rowMeans(
        latent$zero_class == "credible_structural_zero"
      ),
      always_zero_target_fraction = mean(rowSums(rna_counts) == 0),
      zero_pattern_asymmetry = zero_pattern,
      missing_target_fraction = mean(!is.finite(rna_counts)),
      zero_support_sensitive_threshold = 0.25
    ),
    alpha_sensitivity_grid = alpha_grid,
    capacity_params = list(
      regulatory_alpha = regulatory_alpha,
      regulatory_odds_budget = 2^regulatory_alpha,
      gene_half_saturation = gene_half_saturation,
      regulatory_mode = paste0(
        "Pando_outer_OOF_condition_", comparison_support
      ),
      link_function = "tanh(G / shared_oof_scale)",
      promiscuity_mode = "none",
      and_method = gpr_and_method,
      or_method = "sum",
      parallel = parallel
    ),
    projection_provenance = list(
      pando_schema = "pando_condition_grn_fit_v5",
      projection_origin = "outer_condition_stratified_cell_oof",
      projection_used_for_penalty = TRUE,
      full_fit_projection_used_for_penalty = FALSE,
      pando_version = grn_result$pando_installed_version,
      pando_file_fingerprint = grn_result$pando_file_fingerprint,
      supercell_membership = "membership_table(cell_id, metacell_id)",
      projection_order = paste(
        "outer-heldout single-cell TF*ATAC -> training-only balanced",
        "transform -> outer-training coefficient -> SuperCell mean"
      ),
      comparison_support_requested = comparison_support,
      component = "condition",
      unavailable_target_policy = "exclude_from_multiome",
      rna_only_baseline_reported_separately = TRUE
    ),
    inference_class = "metacell_statistical_unit_within_dataset",
    statistical_unit = "metacell",
    metacell_statistical_inference = TRUE,
    biological_replicate_inference = FALSE,
    formal_biological_replicate_pvalue = FALSE,
    pvalue_interpretation =
      "within_dataset_metacell_units_not_sample_level_replicates"
  )
}
