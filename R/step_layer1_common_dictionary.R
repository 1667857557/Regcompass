# Common-dictionary Pando projection used by Layer 1.

.rc_require_pando_condition_grn_fit <- function(fit) {
  if (!inherits(fit, "ConditionGRNFit") ||
      !identical(fit$schema_version, .RC_PANDO_CONDITION_GRN_FIT_SCHEMA) ||
      !identical(fit$projection_effect_column, "penalty_effect") ||
      !identical(fit$projection_policy, "padj_significant_effects_only") ||
      !identical(fit$scale, FALSE) || !identical(fit$interaction, ":")) {
    stop("Pando common-dictionary condition fit is incompatible.",
         call. = FALSE)
  }
  coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
  required <- c(
    "target", "condition", "estimate", "padj", "significant",
    "penalty_effect", "estimable"
  )
  if (!all(required %in% colnames(coefficient)) ||
      any(coefficient$significant & !coefficient$estimable) ||
      any(coefficient$significant &
          (!is.finite(coefficient$padj) |
           coefficient$padj >= as.numeric(fit$padj_threshold))) ||
      any(coefficient$significant &
          coefficient$penalty_effect != coefficient$estimate, na.rm = TRUE) ||
      any(!coefficient$significant & coefficient$penalty_effect != 0,
          na.rm = TRUE)) {
    stop("Pando penalty effects do not match the BH significance contract.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.rc_condition_pando_projection <- function(
    grn_result, membership, unit_meta, genes, comparison_support) {
  primary <- matrix(
    NA_real_, length(genes), nrow(unit_meta),
    dimnames = list(genes, unit_meta$unit_id)
  )
  reliability <- primary
  coverage <- list()
  for (fit in grn_result$condition_grn_fits) {
    .rc_require_pando_condition_grn_fit(fit)
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
    targets <- intersect(rownames(score), rownames(primary))
    units <- intersect(colnames(score), colnames(primary))
    primary[targets, units] <- score[targets, units, drop = FALSE]

    fit_table <- as.data.frame(fit$fit, stringsAsFactors = FALSE)
    fit_table$target <- tolower(as.character(fit_table$target))
    fit_table$condition <- as.character(fit_table$condition)
    condition_col <- as.character(fit$condition_col)
    celltype_col <- as.character(fit$cell_type_col)
    if (!all(c(condition_col, celltype_col) %in% colnames(unit_meta))) {
      stop("Metacell metadata lack fitted condition or cell-type columns.",
           call. = FALSE)
    }
    for (condition in fit$condition_levels) {
      selected_units <- unit_meta$unit_id[
        as.character(unit_meta[[condition_col]]) == condition &
          as.character(unit_meta[[celltype_col]]) == fit$cell_type
      ]
      one <- fit_table[fit_table$condition == condition, , drop = FALSE]
      target <- intersect(one$target, rownames(reliability))
      index <- match(target, one$target)
      value <- sqrt(pmin(1, pmax(0, as.numeric(one$rsq[index]))))
      value[!is.finite(value)] <- NA_real_
      if (length(target) && length(selected_units)) {
        reliability[target, selected_units] <- matrix(
          value, nrow = length(target), ncol = length(selected_units)
        )
      }
    }
    coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
    coverage[[length(coverage) + 1L]] <- data.frame(
      cell_type = fit$cell_type,
      condition = fit$condition_levels,
      n_dictionary_edges = nrow(fit$edge_dictionary),
      n_significant_edges = vapply(fit$condition_levels, function(condition) {
        sum(coefficient$condition == condition &
              coefficient$significant %in% TRUE, na.rm = TRUE)
      }, integer(1)),
      padj_threshold = fit$padj_threshold,
      projection_effect = "penalty_effect",
      stringsAsFactors = FALSE
    )
  }
  common <- primary
  condition_unique <- primary * 0
  list(
    common = common,
    primary = primary,
    condition_unique = condition_unique,
    reliability = reliability,
    coverage = .rc_bind_frames_fill(coverage),
    origin = "paired_cell_full_fit_fixed_dictionary_glm_padj_filtered",
    full_fit_used = TRUE,
    pando_schema = .RC_PANDO_CONDITION_GRN_FIT_SCHEMA,
    primary_projection = "padj_filtered_fixed_dictionary_condition_glm",
    common_projection_role =
      "compatibility_alias_of_primary_no_common_support_decomposition",
    nonestimable_projection_policy =
      "coefficient_NA_and_zero_realized_penalty_contribution",
    comparison_support_ignored = comparison_support
  )
}
