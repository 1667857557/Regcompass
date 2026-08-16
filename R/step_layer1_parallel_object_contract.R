# Route Layer 1 projection to the exact per-cell-type Pando object.

.rc_validate_pando_fit_metadata_frame <- function(
    metadata, fits, condition_col, celltype_col) {
  if (!is.data.frame(metadata) ||
      !all(c(condition_col, celltype_col) %in% colnames(metadata)) ||
      is.null(rownames(metadata)) || anyNA(rownames(metadata)) ||
      any(!nzchar(rownames(metadata))) || anyDuplicated(rownames(metadata))) {
    stop("Pando object metadata cannot validate condition fit cell mappings.",
         call. = FALSE)
  }
  if (inherits(fits, "ConditionGRNFit")) fits <- list(fits)
  if (!is.list(fits) || !length(fits)) {
    stop("No Pando condition fits are available for metadata validation.",
         call. = FALSE)
  }
  for (fit in fits) {
    if (!identical(as.character(fit$condition_col), condition_col) ||
        !identical(as.character(fit$cell_type_col), celltype_col)) {
      stop(
        "Pando fit metadata columns do not match the RegCompass request: ",
        "fit condition_col='", as.character(fit$condition_col),
        "', cell_type_col='", as.character(fit$cell_type_col),
        "'; requested condition_col='", condition_col,
        "', cell_type_col='", celltype_col, "'.", call. = FALSE
      )
    }
    levels <- as.character(fit$condition_levels)
    cells_by_condition <- fit$condition_cell_ids[levels]
    for (condition in levels) {
      cells <- as.character(cells_by_condition[[condition]])
      missing <- setdiff(cells, rownames(metadata))
      if (length(missing)) {
        stop(
          "Pando fit references cells absent from its stored object; first ",
          "missing ID: ", missing[[1L]], ".", call. = FALSE
        )
      }
      observed_condition <- as.character(metadata[cells, condition_col])
      observed_celltype <- as.character(metadata[cells, celltype_col])
      if (anyNA(observed_condition) || anyNA(observed_celltype) ||
          any(observed_condition != condition) ||
          any(observed_celltype != as.character(fit$cell_type))) {
        stop(
          "Pando fit cell assignments disagree with stored object metadata ",
          "for cell type '", as.character(fit$cell_type),
          "' and condition '", condition, "'.", call. = FALSE
        )
      }
    }
  }
  invisible(TRUE)
}

.rc_require_layer1_condition_grn_fit <- function(fit) {
  .rc_require_pando_condition_grn_fit(fit)
  if (is.null(fit$regcompass_penalty_filter)) return(invisible(TRUE))

  threshold <- .rc_condition_padj_threshold(fit = fit)
  rsq_threshold <- .rc_target_rsq_threshold(
    fit$regcompass_target_rsq_threshold %||% .rc_target_rsq_threshold()
  )
  expected_filter <- paste0(
    "Pando BH-active edge & target fit_status == 'ok' & full-data target R2 >= ",
    format(rsq_threshold, trim = TRUE)
  )
  filter_value <- as.character(fit$regcompass_penalty_filter)
  if (length(filter_value) != 1L || is.na(filter_value) ||
      !identical(filter_value, expected_filter)) {
    stop("RegCompass condition-GRN penalty gate metadata are inconsistent.",
         call. = FALSE)
  }

  coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
  required <- c(
    "statistically_supported", "global_support", "local_support",
    "active", "significant", "penalty_effect", "estimate", "estimable",
    "padj", "fit_status", "target_rsq", "target_model_supported",
    "penalty_eligible"
  )
  if (!all(required %in% colnames(coefficient))) {
    stop(
      "RegCompass-gated condition fits require Pando active-edge provenance, ",
      "ridge coefficients, target fit status/R2 and penalty eligibility.",
      call. = FALSE
    )
  }
  .rc_validate_pando_active_condition_edges(
    coefficient, padj_threshold = threshold
  )
  expected_supported <-
    trimws(as.character(coefficient$fit_status)) == "ok" &
    is.finite(as.numeric(coefficient$target_rsq)) &
    as.numeric(coefficient$target_rsq) >= rsq_threshold
  if (!identical(
      as.logical(coefficient$target_model_supported), expected_supported
  )) {
    stop("RegCompass target-model support flags do not match full-data R2 gate.",
         call. = FALSE)
  }
  expected_gate <- .rc_condition_penalty_gate(
    coefficient,
    padj_threshold = threshold,
    target_rsq_threshold = rsq_threshold
  )
  if (!identical(as.logical(coefficient$penalty_eligible), expected_gate)) {
    stop(
      "RegCompass penalty_eligible flags must equal Pando BH-active edges with ",
      "valid target fit status and full-data target R2 support.", call. = FALSE
    )
  }
  invisible(TRUE)
}

.rc_condition_pando_object_for_fit <- function(grn_result, fit) {
  .rc_require_layer1_condition_grn_fit(fit)
  cell_type <- as.character(fit$cell_type)[[1L]]
  object_map <- grn_result$pando_grn_data_by_cell_type %||% list()
  if (length(object_map)) {
    map_names <- names(object_map)
    fit_types <- vapply(grn_result$condition_grn_fits, function(value) {
      as.character(value$cell_type)[[1L]]
    }, character(1))
    if (is.null(map_names) || anyNA(map_names) ||
        any(!nzchar(map_names)) || anyDuplicated(map_names) ||
        anyNA(fit_types) || any(!nzchar(fit_types)) ||
        anyDuplicated(fit_types) || !setequal(map_names, fit_types)) {
      stop(
        "Condition fits and per-cell-type Pando objects do not form the same ",
        "unique cell-type set.", call. = FALSE
      )
    }
    if (!cell_type %in% map_names) {
      stop("No Pando GRNData object is stored for condition-GRN cell type `",
           cell_type, "`.", call. = FALSE)
    }
    object <- object_map[[cell_type]]
  } else {
    object <- grn_result$pando_grn_data
  }
  if (!inherits(object, "GRNData")) {
    stop("Condition-GRN projection requires a Pando GRNData object for `",
         cell_type, "`.", call. = FALSE)
  }
  fit_cells <- as.character(unlist(
    fit$condition_cell_ids[fit$condition_levels], use.names = FALSE
  ))
  object_cells <- colnames(object@data)
  missing <- setdiff(fit_cells, object_cells)
  if (length(missing)) {
    stop(
      "The Pando object for `", cell_type, "` is missing ",
      length(missing), " fitted cell(s); first missing ID: ", missing[[1L]],
      call. = FALSE
    )
  }
  .rc_validate_pando_fit_metadata_frame(
    metadata = object@data@meta.data,
    fits = list(fit),
    condition_col = as.character(fit$condition_col)[[1L]],
    celltype_col = as.character(fit$cell_type_col)[[1L]]
  )
  object
}

.rc_condition_target_reliability <- function(fit) {
  fit_table <- as.data.frame(fit$fit, stringsAsFactors = FALSE)
  coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
  required_fit <- c("target", "condition", "rsq", "fit_status")
  required_coefficient <- c(
    "target", "condition", "active", "estimable", "estimate", "padj",
    "global_support", "local_support"
  )
  if (!nrow(fit_table) ||
      !all(required_fit %in% colnames(fit_table)) ||
      !all(required_coefficient %in% colnames(coefficient))) {
    stop(
      "Condition-GRN target eligibility requires target-condition fit diagnostics ",
      "and Pando active ridge-edge flags.", call. = FALSE
    )
  }

  fit_target <- toupper(trimws(as.character(fit_table$target)))
  fit_condition <- trimws(as.character(fit_table$condition))
  fit_key <- paste(fit_target, fit_condition, sep = "\001")
  if (anyNA(fit_key) || any(!nzchar(fit_target)) ||
      any(!nzchar(fit_condition)) || anyDuplicated(fit_key)) {
    stop(
      "Condition-GRN target-condition fit diagnostics must be complete and unique.",
      call. = FALSE
    )
  }

  coefficient_target <- toupper(trimws(as.character(coefficient$target)))
  coefficient_condition <- trimws(as.character(coefficient$condition))
  coefficient_key <- paste(
    coefficient_target, coefficient_condition, sep = "\001"
  )
  if (anyNA(coefficient_key) || any(!nzchar(coefficient_target)) ||
      any(!nzchar(coefficient_condition))) {
    stop("Condition-GRN coefficients contain incomplete target-condition labels.",
         call. = FALSE)
  }
  coefficient_index <- match(coefficient_key, fit_key)
  if (anyNA(coefficient_index)) {
    stop(
      "Condition-GRN coefficients cannot be aligned to target fit diagnostics.",
      call. = FALSE
    )
  }

  threshold <- .rc_condition_padj_threshold(fit = fit)
  rsq_threshold <- .rc_target_rsq_threshold(
    fit$regcompass_target_rsq_threshold %||% .rc_target_rsq_threshold()
  )
  rsq <- suppressWarnings(as.numeric(fit_table$rsq))
  fit_status <- trimws(as.character(fit_table$fit_status))
  if (anyNA(fit_status) || any(!nzchar(fit_status))) {
    stop("Condition-GRN fit_status values must be complete.", call. = FALSE)
  }

  coefficient$fit_status <- fit_status[coefficient_index]
  coefficient$rsq <- rsq[coefficient_index]
  coefficient$padj_threshold <- threshold
  final_gate <- .rc_condition_penalty_gate(
    coefficient,
    padj_threshold = threshold,
    target_rsq_threshold = rsq_threshold
  )
  n_projection_edges <- tabulate(
    coefficient_index[final_gate], nbins = nrow(fit_table)
  )
  reliability <- rep(NA_real_, nrow(fit_table))
  evaluated <- is.finite(rsq)
  reliability[evaluated] <- as.numeric(
    fit_status[evaluated] == "ok" &
      rsq[evaluated] >= rsq_threshold &
      n_projection_edges[evaluated] > 0L
  )

  data.frame(
    target = as.character(fit_table$target),
    condition = fit_condition,
    rsq = rsq,
    fit_status = fit_status,
    n_significant_edges = as.integer(n_projection_edges),
    n_projection_edges = as.integer(n_projection_edges),
    padj_threshold = threshold,
    target_rsq_threshold = rsq_threshold,
    reliability = reliability,
    stringsAsFactors = FALSE
  )
}

.rc_condition_pando_projection <- function(
    grn_result, membership, unit_meta, genes, rna_assay, atac_assay) {
  projection <- matrix(
    NA_real_, length(genes), nrow(unit_meta),
    dimnames = list(genes, unit_meta$unit_id)
  )
  reliability <- projection
  coverage <- list()

  for (fit_raw in grn_result$condition_grn_fits) {
    fit <- .rc_apply_condition_penalty_gate(fit_raw)
    .rc_require_layer1_condition_grn_fit(fit)
    threshold <- .rc_condition_padj_threshold(fit = fit)
    pando_object <- .rc_condition_pando_object_for_fit(grn_result, fit)

    rna <- Matrix::t(Pando::LayerData(
      pando_object, assay = rna_assay, layer = fit$rna_layer
    ))
    atac <- Matrix::t(Pando::LayerData(
      pando_object, assay = atac_assay, layer = fit$peak_layer
    ))
    cells <- colnames(pando_object@data)
    rownames(rna) <- cells
    rownames(atac) <- cells

    coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
    dictionary <- as.data.frame(fit$edge_dictionary, stringsAsFactors = FALSE)
    coefficient <- merge(
      coefficient,
      dictionary[, c(
        "edge_id", "tf", "target", "region", "atac_feature_id"
      ), drop = FALSE],
      by = "edge_id", all.x = TRUE, sort = FALSE,
      suffixes = c("", ".dictionary")
    )
    projection_gate <- .rc_condition_penalty_gate(
      coefficient,
      padj_threshold = threshold,
      target_rsq_threshold = fit$regcompass_target_rsq_threshold
    )
    coefficient <- coefficient[projection_gate, , drop = FALSE]
    coefficient$estimate <- coefficient$penalty_effect

    fit_units <- unit_meta$unit_id[
      as.character(unit_meta[[fit$cell_type_col]]) == fit$cell_type
    ]
    fit_membership <- membership[
      membership$metacell_id %in% fit_units &
        membership$cell_id %in% cells,
      , drop = FALSE
    ]
    expected_cells <- unique(unlist(
      fit$condition_cell_ids[fit$condition_levels], use.names = FALSE
    ))
    if (anyDuplicated(fit_membership$cell_id) ||
        !setequal(fit_membership$cell_id, expected_cells)) {
      stop(
        "Condition-Pando projection requires exactly one metacell membership ",
        "for every fitted cell.", call. = FALSE
      )
    }

    cells_by_group <- split(
      as.character(fit_membership$cell_id),
      as.character(fit_membership$metacell_id)
    )
    unit_condition <- stats::setNames(
      as.character(unit_meta[[fit$condition_col]]), unit_meta$unit_id
    )
    edges_by_group <- lapply(names(cells_by_group), function(unit) {
      coefficient[
        as.character(coefficient$condition) == unit_condition[[unit]],
        , drop = FALSE
      ]
    })
    names(edges_by_group) <- names(cells_by_group)

    # Canonical RegCompass estimand keeps the paired interaction dominant while
    # shrinking 25% toward the product of separate metacell modality means:
    # beta * [0.75 * mean(TF * ATAC) + 0.25 * mean(TF) * mean(ATAC)].
    score <- .rc_pando_projection_from_group_means(
      rna, atac, edges_by_group, cells_by_group, genes
    )
    targets <- intersect(colnames(score), rownames(projection))
    units <- intersect(rownames(score), colnames(projection))
    if (length(targets) && length(units)) {
      projection[targets, units] <- t(score[units, targets, drop = FALSE])
    }

    target_reliability <- .rc_condition_target_reliability(fit)
    target_reliability$target <- tolower(
      trimws(as.character(target_reliability$target))
    )
    condition_col <- as.character(fit$condition_col)[[1L]]
    celltype_col <- as.character(fit$cell_type_col)[[1L]]
    if (!all(c(condition_col, celltype_col) %in% colnames(unit_meta))) {
      stop("Metacell metadata lack fitted condition or cell-type columns.",
           call. = FALSE)
    }
    for (condition in fit$condition_levels) {
      selected_units <- intersect(
        as.character(unit_meta$unit_id[
          as.character(unit_meta[[condition_col]]) == condition &
            as.character(unit_meta[[celltype_col]]) == fit$cell_type
        ]),
        colnames(reliability)
      )
      one <- target_reliability[
        as.character(target_reliability$condition) == condition &
          is.finite(target_reliability$reliability),
        , drop = FALSE
      ]
      selected_targets <- intersect(one$target, rownames(reliability))
      if (length(selected_targets) && length(selected_units)) {
        q <- one$reliability[match(selected_targets, one$target)]
        reliability[selected_targets, selected_units] <- matrix(
          q,
          nrow = length(selected_targets),
          ncol = length(selected_units)
        )
      }
    }

    all_coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
    coverage[[length(coverage) + 1L]] <- data.frame(
      cell_type = fit$cell_type,
      condition = fit$condition_levels,
      n_dictionary_edges = nrow(fit$edge_dictionary),
      n_projection_edges = vapply(fit$condition_levels, function(condition) {
        sum(coefficient$condition == condition, na.rm = TRUE)
      }, integer(1)),
      n_significant_edges = vapply(fit$condition_levels, function(condition) {
        sum(
          all_coefficient$condition == condition &
            all_coefficient$active %in% TRUE,
          na.rm = TRUE
        )
      }, integer(1)),
      mean_target_reliability = vapply(
        fit$condition_levels,
        function(condition) {
          value <- target_reliability$reliability[
            target_reliability$condition == condition
          ]
          value <- value[is.finite(value)]
          if (length(value)) mean(value) else NA_real_
        },
        numeric(1)
      ),
      reliability_definition = paste0(
        "binary target eligibility after Pando BH edge, fit_status and ",
        "selected-lambda full-data R2 >= ",
        format(fit$regcompass_target_rsq_threshold, trim = TRUE)
      ),
      padj_threshold = threshold,
      target_rsq_threshold = fit$regcompass_target_rsq_threshold,
      corr_threshold = .RC_PANDO_PENALTY_CORR_THRESHOLD,
      estimate_threshold = .RC_PANDO_PENALTY_ESTIMATE_THRESHOLD,
      projection_effect =
        "Pando_condition_BH_active_no_fusion_ridge_penalty_effect",
      pando_object_scope = "cell_type_exact_feature_space",
      aggregation_contract =
        "beta_times_0.75_paired_mean_product_plus_0.25_product_of_means",
      paired_interaction_weight =
        1 - .RC_PANDO_PROJECTION_PRODUCT_OF_MEANS_WEIGHT,
      product_of_means_weight =
        .RC_PANDO_PROJECTION_PRODUCT_OF_MEANS_WEIGHT,
      stringsAsFactors = FALSE
    )
  }

  list(
    projection = projection,
    reliability = reliability,
    coverage = .rc_bind_frames_fill(coverage),
    origin = "paired_cell_common_dictionary_condition_bh_active_edges",
    pando_schema = .RC_PANDO_CONDITION_GRN_FIT_SCHEMA,
    projection_name = "active_no_fusion_ridge_condition_effect",
    nonestimable_policy =
      "nonactive_or_nonestimable_condition_edge_has_zero_projection_contribution"
  )
}
