# Canonical sign-only Pando projection used by Layer 1.
#
# Pando coefficient magnitudes remain available for inference diagnostics and
# auditing, but they are not used as RegCompass penalty magnitudes.  The fitted
# coefficient contributes only its sign after the existing estimability,
# target-fit-status and BH gates.  To avoid making the sign-only score depend on
# arbitrary RNA/ATAC feature units, each paired-cell TF-by-ATAC predictor is
# calibrated by a positive, cell-type-wide robust scale before membership
# aggregation.

.rc_sign_only_edge_direction <- function(estimate, eligible = NULL) {
  estimate <- suppressWarnings(as.numeric(estimate))
  if (is.null(eligible)) eligible <- rep(TRUE, length(estimate))
  eligible <- as.logical(eligible)
  if (length(eligible) != length(estimate) || anyNA(eligible)) {
    stop("Sign-only Pando edge eligibility must align to coefficient estimates.",
         call. = FALSE)
  }
  answer <- numeric(length(estimate))
  keep <- eligible & is.finite(estimate) & estimate != 0
  answer[keep] <- sign(estimate[keep])
  answer
}

.rc_sign_only_predictor_activity <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!length(x) || any(!is.finite(x))) {
    stop("Sign-only Pando predictor activity requires finite paired-cell values.",
         call. = FALSE)
  }
  tolerance <- sqrt(.Machine$double.eps) * pmax(1, abs(x))
  if (any(x < -tolerance)) {
    stop(
      "Sign-only Pando predictor activity requires non-negative normalized ",
      "RNA-by-ATAC products.", call. = FALSE
    )
  }
  x <- pmax(x, 0)
  scale <- if (length(x) >= 2L && any(x > 0)) {
    max(
      stats::IQR(x) / 1.349,
      stats::mad(x, constant = 1.4826),
      sqrt(mean(x^2)),
      1e-6,
      na.rm = TRUE
    )
  } else {
    1
  }
  if (!is.finite(scale) || scale <= 0) scale <- 1
  list(
    activity = tanh(x / scale),
    scale = scale
  )
}

.rc_sign_only_pair_key <- function(tf, region) {
  paste(
    toupper(trimws(as.character(tf))),
    .rc_pando_region_key(as.character(region)),
    sep = "\001"
  )
}

.rc_sign_only_pair_metacell_activity <- function(
    rna, atac, edge, membership, cells) {
  if (is.null(dim(rna)) || is.null(dim(atac)) ||
      is.null(rownames(rna)) || is.null(rownames(atac)) ||
      is.null(colnames(rna)) || is.null(colnames(atac))) {
    stop("Sign-only Pando projection requires named RNA and ATAC matrices.",
         call. = FALSE)
  }
  edge <- as.data.frame(edge, stringsAsFactors = FALSE)
  if (!all(c("tf", "region") %in% colnames(edge))) {
    stop("Sign-only Pando projection requires TF and region edge columns.",
         call. = FALSE)
  }
  if (!is.data.frame(membership) ||
      !all(c("cell_id", "metacell_id") %in% colnames(membership)) ||
      anyDuplicated(as.character(membership$cell_id))) {
    stop("Sign-only Pando projection requires unique exact metacell membership.",
         call. = FALSE)
  }
  cells <- as.character(cells)
  if (!length(cells) || anyNA(cells) || any(!nzchar(cells)) ||
      anyDuplicated(cells)) {
    stop("Sign-only Pando projection cells must be unique and non-empty.",
         call. = FALSE)
  }
  if (anyNA(match(cells, rownames(rna))) || anyNA(match(cells, rownames(atac)))) {
    stop("Pando RNA/ATAC matrices do not contain every fitted projection cell.",
         call. = FALSE)
  }
  membership_index <- match(cells, as.character(membership$cell_id))
  if (anyNA(membership_index)) {
    stop("Every fitted Pando cell must have exact SuperCell membership.",
         call. = FALSE)
  }
  group <- as.character(membership$metacell_id[membership_index])
  if (anyNA(group) || any(!nzchar(group))) {
    stop("Metacell membership IDs must be complete for sign-only projection.",
         call. = FALSE)
  }
  units <- unique(group)
  group_index <- match(group, units)
  group_n <- tabulate(group_index, nbins = length(units))

  pair_key <- .rc_sign_only_pair_key(edge$tf, edge$region)
  pair_rows <- !duplicated(pair_key)
  pair <- edge[pair_rows, c("tf", "region"), drop = FALSE]
  pair$pair_key <- pair_key[pair_rows]
  pair$tf_key <- toupper(trimws(as.character(pair$tf)))
  pair$region_key <- .rc_pando_region_key(as.character(pair$region))

  rna_key <- toupper(colnames(rna))
  atac_key <- .rc_pando_region_key(colnames(atac))
  if (anyDuplicated(rna_key)) {
    stop("RNA features are duplicated after case normalization.", call. = FALSE)
  }
  if (anyDuplicated(atac_key)) {
    stop("ATAC features are duplicated after Pando region normalization.",
         call. = FALSE)
  }
  tf_index <- match(pair$tf_key, rna_key)
  peak_index <- match(pair$region_key, atac_key)
  usable <- !is.na(tf_index) & !is.na(peak_index)
  pair <- pair[usable, , drop = FALSE]
  tf_index <- tf_index[usable]
  peak_index <- peak_index[usable]

  activity <- matrix(
    0,
    nrow = nrow(pair), ncol = length(units),
    dimnames = list(as.character(pair$pair_key), units)
  )
  scale <- numeric(nrow(pair))
  if (nrow(pair)) {
    rna_cell_index <- match(cells, rownames(rna))
    atac_cell_index <- match(cells, rownames(atac))
    group_factor <- factor(group_index, levels = seq_along(units))
    for (i in seq_len(nrow(pair))) {
      predictor <-
        as.numeric(rna[rna_cell_index, tf_index[[i]]]) *
        as.numeric(atac[atac_cell_index, peak_index[[i]]])
      transformed <- .rc_sign_only_predictor_activity(predictor)
      grouped_sum <- rowsum(
        matrix(transformed$activity, ncol = 1L),
        group = group_factor,
        reorder = FALSE
      )
      activity[i, ] <- as.numeric(grouped_sum[, 1L]) / group_n
      scale[[i]] <- transformed$scale
    }
  }

  list(
    activity = activity,
    scale = data.frame(
      pair_key = as.character(pair$pair_key),
      tf = as.character(pair$tf),
      region = as.character(pair$region),
      predictor_scale = as.numeric(scale),
      stringsAsFactors = FALSE
    ),
    n_requested_pairs = length(unique(pair_key)),
    n_mapped_pairs = nrow(pair),
    units = units
  )
}

.rc_condition_pando_projection_sign_only <- function(
    grn_result, membership, unit_meta, genes) {
  projection <- matrix(
    NA_real_, length(genes), nrow(unit_meta),
    dimnames = list(genes, unit_meta$unit_id)
  )
  reliability <- projection
  coverage <- list()

  for (fit in grn_result$condition_grn_fits) {
    .rc_require_layer1_condition_grn_fit(fit)
    pando_object <- .rc_condition_pando_object_for_fit(grn_result, fit)
    cells <- as.character(unlist(
      fit$condition_cell_ids[fit$condition_levels], use.names = FALSE
    ))
    pando_params <- Pando::Params(pando_object)
    condition_rna_assay <- as.character(pando_params$rna_assay %||% "RNA")[[1L]]
    condition_atac_assay <- as.character(pando_params$peak_assay %||% "ATAC")[[1L]]
    rna <- Matrix::t(Pando::LayerData(
      pando_object, assay = condition_rna_assay, layer = "data"
    ))
    atac <- Matrix::t(Pando::LayerData(
      pando_object, assay = condition_atac_assay, layer = "data"
    ))
    rownames(rna) <- colnames(pando_object@data)
    rownames(atac) <- colnames(pando_object@data)

    coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
    gate <- .rc_condition_penalty_gate(coefficient)
    coefficient$penalty_direction <- .rc_sign_only_edge_direction(
      coefficient$estimate, gate
    )
    coefficient$target_key <- tolower(trimws(as.character(coefficient$target)))
    coefficient$pair_key <- .rc_sign_only_pair_key(
      coefficient$tf, coefficient$region
    )
    active <- coefficient[
      gate & coefficient$penalty_direction != 0,
      , drop = FALSE
    ]

    condition_col <- as.character(fit$condition_col)[[1L]]
    celltype_col <- as.character(fit$cell_type_col)[[1L]]
    if (!all(c(condition_col, celltype_col) %in% colnames(unit_meta))) {
      stop("Metacell metadata lack fitted condition or cell-type columns.",
           call. = FALSE)
    }
    fit_units <- unit_meta$unit_id[
      as.character(unit_meta[[celltype_col]]) == as.character(fit$cell_type)
    ]
    fit_targets <- intersect(
      tolower(as.character(fit$target_genes)), rownames(projection)
    )
    if (length(fit_targets) && length(fit_units)) {
      projection[fit_targets, fit_units] <- 0
    }

    pair_projection <- if (nrow(active)) {
      .rc_sign_only_pair_metacell_activity(
        rna = rna,
        atac = atac,
        edge = active,
        membership = membership,
        cells = cells
      )
    } else {
      list(
        activity = matrix(0, 0, 0),
        scale = data.frame(),
        n_requested_pairs = 0L,
        n_mapped_pairs = 0L,
        units = character()
      )
    }

    for (condition in fit$condition_levels) {
      selected_units <- unit_meta$unit_id[
        as.character(unit_meta[[condition_col]]) == condition &
          as.character(unit_meta[[celltype_col]]) == as.character(fit$cell_type)
      ]
      one <- active[
        as.character(active$condition) == condition,
        , drop = FALSE
      ]
      mapped_rows <- one$pair_key %in% rownames(pair_projection$activity) &
        one$target_key %in% rownames(projection)
      one <- one[mapped_rows, , drop = FALSE]
      mapped_units <- intersect(
        selected_units, colnames(pair_projection$activity)
      )
      if (nrow(one) && length(mapped_units)) {
        for (i in seq_len(nrow(one))) {
          target <- one$target_key[[i]]
          pair <- one$pair_key[[i]]
          projection[target, mapped_units] <-
            projection[target, mapped_units] +
            as.numeric(one$penalty_direction[[i]]) *
            pair_projection$activity[pair, mapped_units]
        }
        reliable_targets <- intersect(
          unique(as.character(one$target_key)), rownames(reliability)
        )
        reliability[reliable_targets, mapped_units] <- 1
      }
    }

    coverage[[length(coverage) + 1L]] <- data.frame(
      cell_type = as.character(fit$cell_type),
      condition = as.character(fit$condition_levels),
      n_dictionary_edges = nrow(fit$edge_dictionary),
      n_significant_edges = vapply(fit$condition_levels, function(condition) {
        sum(
          gate & as.character(coefficient$condition) == condition,
          na.rm = TRUE
        )
      }, integer(1)),
      n_signed_penalty_edges = vapply(fit$condition_levels, function(condition) {
        sum(
          coefficient$penalty_direction != 0 &
            as.character(coefficient$condition) == condition,
          na.rm = TRUE
        )
      }, integer(1)),
      n_requested_predictor_pairs = pair_projection$n_requested_pairs,
      n_mapped_predictor_pairs = pair_projection$n_mapped_pairs,
      padj_threshold = 0.05,
      corr_threshold = .RC_PANDO_PENALTY_CORR_THRESHOLD,
      estimate_threshold = .RC_PANDO_PENALTY_ESTIMATE_THRESHOLD,
      estimate_magnitude_used_for_penalty = FALSE,
      projection_effect =
        "sign(estimate)*tanh((paired_cell_TF*ATAC)/celltype_edge_scale)",
      pando_object_scope = "cell_type_exact_feature_space",
      aggregation_contract = "all_projected_cells_have_exact_membership",
      stringsAsFactors = FALSE
    )
  }

  list(
    projection = projection,
    reliability = reliability,
    coverage = .rc_bind_frames_fill(coverage),
    origin = "paired_cell_sign_only_fixed_dictionary_glm_bh_filtered",
    pando_schema = .RC_PANDO_CONDITION_GRN_FIT_SCHEMA,
    projection_name = "bh_filtered_sign_only_scaled_pair_activity",
    nonestimable_policy =
      "coefficient_NA_and_zero_realized_penalty_contribution",
    estimate_magnitude_used_for_penalty = FALSE
  )
}

.rc_standard_pando_projection_sign_only <- function(
    grn_result, membership, metacell_meta, condition_col, celltype_col,
    rna_assay, atac_assay, target_genes) {
  units <- as.character(metacell_meta$unit_id %||% metacell_meta$metacell_id)
  projection <- matrix(
    NA_real_, length(target_genes), length(units),
    dimnames = list(tolower(target_genes), units)
  )
  reliability <- projection
  coverage <- list()

  for (celltype in names(grn_result$standard_pando_objects)) {
    grn <- grn_result$standard_pando_objects[[celltype]]
    cells <- colnames(grn@data)
    one_membership <- membership[
      membership$cell_id %in% cells, , drop = FALSE
    ]
    if (!nrow(one_membership)) next
    rna <- Matrix::t(Pando::LayerData(
      grn, assay = rna_assay, layer = "data"
    ))
    atac <- Matrix::t(Pando::LayerData(
      grn, assay = atac_assay, layer = "data"
    ))
    rownames(rna) <- colnames(grn@data)
    rownames(atac) <- colnames(grn@data)

    all_edge <- grn_result$tf_peak_gene_condition_all[
      as.character(grn_result$tf_peak_gene_condition_all[[celltype_col]]) ==
        celltype,
      , drop = FALSE
    ]
    edge <- grn_result$tf_peak_gene_condition[
      as.character(grn_result$tf_peak_gene_condition[[celltype_col]]) == celltype,
      , drop = FALSE
    ]
    fit_units <- metacell_meta$unit_id[
      as.character(metacell_meta[[celltype_col]]) == celltype
    ]
    fit_targets <- intersect(
      tolower(unique(as.character(all_edge$target))), rownames(projection)
    )
    if (length(fit_targets) && length(fit_units)) {
      projection[fit_targets, fit_units] <- 0
    }
    if (!nrow(edge)) next

    edge$target_key <- tolower(trimws(as.character(edge$target)))
    edge$pair_key <- .rc_sign_only_pair_key(edge$tf, edge$region)
    edge$penalty_direction <- .rc_sign_only_edge_direction(
      edge$estimate, rep(TRUE, nrow(edge))
    )
    edge <- edge[
      edge$penalty_direction != 0 &
        edge$target_key %in% rownames(projection),
      , drop = FALSE
    ]
    if (!nrow(edge)) next

    pair_projection <- .rc_sign_only_pair_metacell_activity(
      rna = rna,
      atac = atac,
      edge = edge,
      membership = one_membership,
      cells = cells
    )
    mapped_units <- intersect(fit_units, colnames(pair_projection$activity))
    mapped_rows <- edge$pair_key %in% rownames(pair_projection$activity)
    edge <- edge[mapped_rows, , drop = FALSE]
    if (!nrow(edge) || !length(mapped_units)) next

    for (i in seq_len(nrow(edge))) {
      target <- edge$target_key[[i]]
      pair <- edge$pair_key[[i]]
      projection[target, mapped_units] <-
        projection[target, mapped_units] +
        as.numeric(edge$penalty_direction[[i]]) *
        pair_projection$activity[pair, mapped_units]
    }

    target <- intersect(unique(edge$target_key), rownames(projection))
    rsq <- tapply(as.numeric(edge$rsq), edge$target_key, function(x) {
      x <- x[is.finite(x)]
      if (length(x)) max(x) else NA_real_
    })
    q <- sqrt(pmin(1, pmax(0, rsq[target])))
    reliability[target, mapped_units] <- matrix(
      q, nrow = length(target), ncol = length(mapped_units)
    )
    coverage[[celltype]] <- data.frame(
      target = target,
      cell_type = celltype,
      n_standard_pando_edges = vapply(target, function(gene) {
        sum(edge$target_key == gene)
      }, integer(1)),
      n_requested_predictor_pairs = pair_projection$n_requested_pairs,
      n_mapped_predictor_pairs = pair_projection$n_mapped_pairs,
      projection_origin = "standard_pando_sign_only_scaled_pair_activity",
      projection_used_for_penalty = TRUE,
      estimate_magnitude_used_for_penalty = FALSE,
      condition_coefficients_calculated = FALSE,
      stringsAsFactors = FALSE
    )
  }

  list(
    projection = projection,
    reliability = reliability,
    coverage = .rc_bind_frames_fill(coverage),
    projection_origin = "standard_pando_sign_only_scaled_pair_activity",
    projection_name = "standard_pando_sign_only_scaled_pair_activity",
    projection_used_for_penalty = TRUE,
    full_fit_projection_used_for_penalty = FALSE,
    estimate_magnitude_used_for_penalty = FALSE
  )
}
