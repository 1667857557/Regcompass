# Condition-full Pando contract and result helpers.

.rc_require_pando_condition_grn_fit <- function(fit) {
  if (!inherits(fit, "ConditionGRNFit") ||
      !identical(fit$schema_version, .RC_PANDO_CONDITION_GRN_FIT_SCHEMA)) {
    stop(
      "RegCompass requires the canonical pando_condition_grn_fit contract.",
      call. = FALSE
    )
  }
  required_masks <- c(
    "coefficient_estimable_mask",
    "projectable_structural_zero_mask",
    "projection_support_mask"
  )
  if (!all(required_masks %in% names(fit))) {
    stop(
      "Pando fit lacks condition-full structural-zero projection masks.",
      call. = FALSE
    )
  }
  reference <- as.matrix(fit$estimability_mask)
  masks <- lapply(required_masks, function(name) as.matrix(fit[[name]]))
  names(masks) <- required_masks
  if (any(vapply(masks, function(value) {
        !is.logical(value) || anyNA(value) ||
          !identical(dimnames(value), dimnames(reference))
      }, logical(1)))) {
    stop("Pando condition-full projection masks are invalid.",
         call. = FALSE)
  }
  expected_zero <- as.matrix(fit$topology_mask) & !reference
  expected_support <- reference | expected_zero
  if (!identical(masks$coefficient_estimable_mask, reference) ||
      !identical(masks$projectable_structural_zero_mask, expected_zero) ||
      !identical(masks$projection_support_mask, expected_support) ||
      !identical(fit$primary_projection, "projection_condition_full_oof") ||
      !identical(
        fit$nonestimable_projection_policy,
        "structural_zero_by_condition"
      )) {
    stop("Pando condition-full projection contract is inconsistent.",
         call. = FALSE)
  }
  full <- as.matrix(fit$projection_condition_full_oof)
  common <- as.matrix(fit$projection_common_oof)
  if (!identical(dimnames(full), dimnames(common)) ||
      any(!is.finite(full)) || any(!is.finite(common))) {
    stop("Pando OOF projections must be aligned and finite.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.rc_extract_condition_grn_contract_unmodified <-
  .rc_extract_condition_grn_contract

.rc_extract_condition_grn_contract <- function(
    grn_object, condition_col, celltype_col,
    min_abs_estimate = 0, min_model_rsq = 0.1) {
  answer <- .rc_extract_condition_grn_contract_unmodified(
    grn_object = grn_object,
    condition_col = condition_col,
    celltype_col = celltype_col,
    min_abs_estimate = min_abs_estimate,
    min_model_rsq = min_model_rsq
  )
  if (!nrow(answer$condition_all)) return(answer)
  answer$condition_all$coefficient_estimable <- FALSE
  answer$condition_all$projectable_structural_zero <- FALSE
  answer$condition_all$projection_supported <- FALSE
  for (fit in answer$fit_contracts) {
    edge_id <- as.character(fit$edge_table$edge_id)
    for (condition in fit$condition_levels) {
      rows <- answer$condition_all$edge_id %in% edge_id &
        as.character(answer$condition_all[[condition_col]]) == condition &
        as.character(answer$condition_all[[celltype_col]]) == fit$cell_type
      index <- match(answer$condition_all$edge_id[rows], edge_id)
      answer$condition_all$coefficient_estimable[rows] <-
        fit$coefficient_estimable_mask[index, condition]
      answer$condition_all$projectable_structural_zero[rows] <-
        fit$projectable_structural_zero_mask[index, condition]
      answer$condition_all$projection_supported[rows] <-
        fit$projection_support_mask[index, condition]
    }
  }
  answer$condition_effect_all <- answer$condition_all
  active_key <- paste(
    answer$condition_active$edge_id,
    answer$condition_active[[condition_col]],
    answer$condition_active[[celltype_col]],
    sep = "\001"
  )
  all_key <- paste(
    answer$condition_all$edge_id,
    answer$condition_all[[condition_col]],
    answer$condition_all[[celltype_col]],
    sep = "\001"
  )
  active_index <- match(active_key, all_key)
  answer$condition_active$coefficient_estimable <-
    answer$condition_all$coefficient_estimable[active_index]
  answer$condition_active$projectable_structural_zero <-
    answer$condition_all$projectable_structural_zero[active_index]
  answer$condition_active$projection_supported <-
    answer$condition_all$projection_supported[active_index]
  answer$condition_effect_active <- answer$condition_active
  answer$projection_contract <- list(
    primary = "condition_full_oof",
    common_component = "common_support_oof",
    condition_unique_component =
      "condition_full_oof_minus_common_support_oof",
    nonestimable_edge = "projectable_structural_zero_by_condition"
  )
  answer
}

.rc_condition_increment_summary <- function(
    microcompass, condition_col, celltype_col) {
  increment <- as.matrix(
    microcompass$penalty_condition_unique_increment
  )
  vmax <- as.matrix(microcompass$vmax)
  if (!identical(dimnames(increment), dimnames(vmax))) {
    stop("Condition-unique increment and vmax matrices are not aligned.",
         call. = FALSE)
  }
  omega <- microcompass$params$omega %||% 0.95
  normalized <- increment / (omega * vmax)
  normalized[
    !is.finite(increment) | !is.finite(vmax) | vmax <= 0
  ] <- NA_real_
  meta <- microcompass$unit_meta
  id_col <- if ("unit_id" %in% colnames(meta)) "unit_id" else "pool_id"
  unit_id <- as.character(meta[[id_col]])
  if (!setequal(unit_id, colnames(normalized))) {
    stop("Condition-unique increment units and metadata differ.",
         call. = FALSE)
  }
  meta <- meta[match(colnames(normalized), unit_id), , drop = FALSE]
  row_meta <- rc_parse_microcompass_row_id(rownames(normalized))
  strata <- unique(meta[, c(condition_col, celltype_col), drop = FALSE])
  rows <- lapply(seq_len(nrow(strata)), function(i) {
    selected <- as.character(meta[[condition_col]]) ==
      as.character(strata[[condition_col]][[i]]) &
      as.character(meta[[celltype_col]]) ==
      as.character(strata[[celltype_col]][[i]])
    value <- matrixStats::rowMedians(
      normalized[, selected, drop = FALSE],
      na.rm = TRUE
    )
    value[is.nan(value)] <- NA_real_
    data.frame(
      row_id = rownames(normalized),
      reaction_id = row_meta$reaction_id,
      target_direction = row_meta$target_direction,
      medium_scenario = row_meta$medium_scenario,
      condition = as.character(strata[[condition_col]][[i]]),
      cell_type = as.character(strata[[celltype_col]][[i]]),
      median_condition_unique_increment_per_required_target_flux = value,
      n_metacells = sum(selected),
      inference_class = "metacell_statistical_unit_within_dataset",
      statistical_unit = "metacell",
      metacell_statistical_inference = TRUE,
      biological_replicate_inference = FALSE,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
