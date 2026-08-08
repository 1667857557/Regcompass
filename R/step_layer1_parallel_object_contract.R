# Route Layer 1 projection to the exact per-cell-type Pando object.

.rc_require_layer1_condition_grn_fit <- function(fit) {
  if (is.null(fit$regcompass_penalty_filter)) {
    .rc_require_pando_condition_grn_fit(fit)
    return(invisible(TRUE))
  }

  expected_filter <-
    "estimable & BH padj < 0.05 & abs(corr) >= 0.05 & abs(estimate) >= 0.05"
  filter_value <- as.character(fit$regcompass_penalty_filter)
  corr_threshold <- suppressWarnings(as.numeric(
    fit$regcompass_corr_threshold
  ))
  estimate_threshold <- suppressWarnings(as.numeric(
    fit$regcompass_estimate_threshold
  ))
  if (length(filter_value) != 1L || is.na(filter_value) ||
      !identical(filter_value, expected_filter) ||
      length(corr_threshold) != 1L ||
      !isTRUE(all.equal(
        corr_threshold, .RC_PANDO_PENALTY_CORR_THRESHOLD
      )) ||
      length(estimate_threshold) != 1L ||
      !isTRUE(all.equal(
        estimate_threshold, .RC_PANDO_PENALTY_ESTIMATE_THRESHOLD
      ))) {
    stop(
      "RegCompass condition-GRN penalty gate metadata are inconsistent.",
      call. = FALSE
    )
  }

  coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
  if (!all(c("significant", "penalty_effect") %in%
           colnames(coefficient))) {
    stop(
      "RegCompass-gated condition fits require significant and penalty_effect.",
      call. = FALSE
    )
  }
  expected_gate <- .rc_condition_penalty_gate(coefficient)
  observed_gate <- as.logical(coefficient$significant)
  if (!identical(observed_gate, expected_gate)) {
    stop(
      "RegCompass-gated significant-edge flags do not match the final penalty gate.",
      call. = FALSE
    )
  }

  estimate <- suppressWarnings(as.numeric(coefficient$estimate))
  expected_effect <- ifelse(expected_gate, estimate, 0)
  observed_effect <- suppressWarnings(as.numeric(
    coefficient$penalty_effect
  ))
  comparable <- is.finite(expected_effect) & is.finite(observed_effect)
  if (any(is.finite(expected_effect) != is.finite(observed_effect)) ||
      any(abs(expected_effect[comparable] -
              observed_effect[comparable]) > 1e-12)) {
    stop(
      "RegCompass-gated penalty_effect does not match the final penalty gate.",
      call. = FALSE
    )
  }

  validation_fit <- fit
  padj <- suppressWarnings(as.numeric(coefficient$padj))
  pando_significant <- coefficient$estimable %in% TRUE &
    is.finite(padj) & padj < 0.05
  coefficient$significant <- pando_significant
  coefficient$penalty_effect <- ifelse(
    pando_significant, estimate, 0
  )
  validation_fit$coefficients <- coefficient
  .rc_require_pando_condition_grn_fit(validation_fit)
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

.rc_condition_pando_projection <- function(
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
    cell_projection <- Pando::project_condition_grn_cells(
      object = pando_object,
      fit = fit,
      targets = genes,
      significant_only = TRUE,
      return_edge_contributions = FALSE
    )
    projection_membership <-
      Pando::validate_condition_grn_projection_membership(
        cell_projection, membership, group_col = "metacell_id"
      )
    aggregated <- Pando::aggregate_condition_grn_projection(
      cell_projection, projection_membership, group_col = "metacell_id"
    )
    score <- as.matrix(aggregated$gene_score)
    rownames(score) <- tolower(rownames(score))
    targets <- intersect(rownames(score), rownames(projection))
    units <- intersect(colnames(score), colnames(projection))
    projection[targets, units] <- score[targets, units, drop = FALSE]

    coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
    condition_col <- as.character(fit$condition_col)[[1L]]
    celltype_col <- as.character(fit$cell_type_col)[[1L]]
    if (!all(c(condition_col, celltype_col) %in% colnames(unit_meta))) {
      stop("Metacell metadata lack fitted condition or cell-type columns.",
           call. = FALSE)
    }
    for (condition in fit$condition_levels) {
      selected_units <- unit_meta$unit_id[
        as.character(unit_meta[[condition_col]]) == condition &
          as.character(unit_meta[[celltype_col]]) == fit$cell_type
      ]
      significant_edges <- coefficient[
        as.character(coefficient$condition) == condition &
          coefficient$significant %in% TRUE,
        , drop = FALSE
      ]
      reliable_targets <- intersect(
        tolower(unique(as.character(significant_edges$target))),
        rownames(reliability)
      )
      if (length(reliable_targets) && length(selected_units)) {
        reliability[reliable_targets, selected_units] <- 1
      }
    }
    coverage[[length(coverage) + 1L]] <- data.frame(
      cell_type = fit$cell_type,
      condition = fit$condition_levels,
      n_dictionary_edges = nrow(fit$edge_dictionary),
      n_significant_edges = vapply(fit$condition_levels, function(condition) {
        sum(
          coefficient$condition == condition &
            coefficient$significant %in% TRUE,
          na.rm = TRUE
        )
      }, integer(1)),
      padj_threshold = 0.05,
      corr_threshold = .RC_PANDO_PENALTY_CORR_THRESHOLD,
      estimate_threshold = .RC_PANDO_PENALTY_ESTIMATE_THRESHOLD,
      projection_effect = "penalty_effect",
      pando_object_scope = "cell_type_exact_feature_space",
      aggregation_contract = "all_projected_cells_have_exact_membership",
      stringsAsFactors = FALSE
    )
  }

  list(
    projection = projection,
    reliability = reliability,
    coverage = .rc_bind_frames_fill(coverage),
    origin = "paired_cell_fixed_dictionary_glm_padj_corr_effect_filtered",
    pando_schema = .RC_PANDO_CONDITION_GRN_FIT_SCHEMA,
    projection_name =
      "padj_corr_effect_filtered_fixed_dictionary_condition_glm",
    nonestimable_policy =
      "coefficient_NA_and_zero_realized_penalty_contribution"
  )
}
