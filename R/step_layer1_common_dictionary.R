# Common-dictionary Pando projection used by Layer 1.

.rc_condition_pando_projection <- function(
    grn_result, membership, unit_meta, genes) {
  projection <- matrix(
    NA_real_, length(genes), nrow(unit_meta),
    dimnames = list(genes, unit_meta$unit_id)
  )
  reliability <- projection
  coverage <- list()

  for (fit in grn_result$condition_grn_fits) {
    fit <- .rc_apply_condition_penalty_gate(fit)
    cell_projection <- Pando::project_condition_grn_cells(
      object = grn_result$pando_grn_data,
      fit = fit,
      targets = genes,
      significant_only = TRUE,
      return_edge_contributions = FALSE
    )
    aggregated <- Pando::aggregate_condition_grn_projection(
      cell_projection, membership, group_col = "metacell_id"
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
