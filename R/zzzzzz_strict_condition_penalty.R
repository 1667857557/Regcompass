# Final common-dictionary contract loaded after the core workflow.
# The multi-condition penalty uses only estimable BH padj < 0.05 coefficients.

.rc_require_pando_condition_grn_fit_strict_base <-
  .rc_require_pando_condition_grn_fit

.rc_require_pando_condition_grn_fit <- function(fit) {
  .rc_require_pando_condition_grn_fit_strict_base(fit)
  if (!identical(toupper(as.character(fit$adjust_method)), "BH") ||
      !isTRUE(all.equal(as.numeric(fit$padj_threshold), 0.05))) {
    stop("RegCompass condition GRNs require BH padj < 0.05.",
         call. = FALSE)
  }
  edge <- as.data.frame(fit$edge_dictionary, stringsAsFactors = FALSE)
  coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
  required <- c(
    "edge_id", "target", "tf", "region", "condition", "estimate",
    "pval", "padj", "significant", "penalty_effect", "estimable"
  )
  if (!all(required %in% colnames(coefficient)) || !nrow(edge) ||
      anyDuplicated(edge$edge_id)) {
    stop("Pando common-dictionary coefficient table is incomplete.",
         call. = FALSE)
  }
  for (condition in fit$condition_levels) {
    one <- coefficient[
      as.character(coefficient$condition) == condition, , drop = FALSE
    ]
    if (!identical(sort(as.character(one$edge_id)),
                   sort(as.character(edge$edge_id)))) {
      stop(
        "Every condition must contain every frozen dictionary edge exactly once.",
        call. = FALSE
      )
    }
    valid <- is.finite(as.numeric(one$pval))
    expected_padj <- rep(NA_real_, nrow(one))
    expected_padj[valid] <- stats::p.adjust(
      as.numeric(one$pval[valid]), method = "BH"
    )
    comparable <- is.finite(expected_padj) & is.finite(as.numeric(one$padj))
    if (any(is.finite(expected_padj) != is.finite(as.numeric(one$padj))) ||
        any(abs(expected_padj[comparable] -
                as.numeric(one$padj[comparable])) > 1e-10)) {
      stop("Stored condition padj values do not equal BH adjustment.",
           call. = FALSE)
    }
  }
  expected_significant <- coefficient$estimable %in% TRUE &
    is.finite(as.numeric(coefficient$padj)) &
    as.numeric(coefficient$padj) < 0.05
  if (!identical(as.logical(coefficient$significant), expected_significant)) {
    stop("Pando significant-edge flags are not exactly estimable & padj < 0.05.",
         call. = FALSE)
  }
  expected_effect <- ifelse(
    expected_significant, as.numeric(coefficient$estimate), 0
  )
  observed_effect <- as.numeric(coefficient$penalty_effect)
  if (any(is.finite(expected_effect) != is.finite(observed_effect)) ||
      any(abs(expected_effect[is.finite(expected_effect)] -
              observed_effect[is.finite(expected_effect)]) > 1e-12)) {
    stop("Pando penalty_effect does not match the strict BH edge contract.",
         call. = FALSE)
  }
  if ("direction" %in% colnames(coefficient)) {
    expected_direction <- ifelse(
      !coefficient$estimable, "undefined",
      ifelse(coefficient$estimate > 0, "positive",
             ifelse(coefficient$estimate < 0, "negative", "zero"))
    )
    if (!identical(as.character(coefficient$direction), expected_direction)) {
      stop("Pando coefficient directions are inconsistent with estimates.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

.rc_extract_condition_grn_contract_strict_base <-
  .rc_extract_condition_grn_contract

.rc_extract_condition_grn_contract <- function(
    grn_object, condition_col, celltype_col,
    min_abs_estimate = 0, min_model_rsq = 0.1) {
  answer <- .rc_extract_condition_grn_contract_strict_base(
    grn_object = grn_object,
    condition_col = condition_col,
    celltype_col = celltype_col,
    min_abs_estimate = 0,
    min_model_rsq = -1e100
  )
  all_edges <- answer$condition_all
  strict <- all_edges$estimable %in% TRUE &
    is.finite(as.numeric(all_edges$padj)) &
    as.numeric(all_edges$padj) < 0.05 &
    is.finite(as.numeric(all_edges$penalty_effect))
  all_edges$penalty_eligible <- strict
  all_edges$active_in_condition <- strict
  answer$condition_all <- all_edges
  answer$condition_effect_all <- all_edges
  answer$condition_active <- all_edges[strict, , drop = FALSE]
  answer$condition_effect_active <- answer$condition_active
  answer$active_tol <- 0
  answer$penalty_filter <-
    "estimable & BH-adjusted padj < 0.05; no effect-size or model-R2 gate"
  answer
}

.rc_condition_pando_projection_strict_base <-
  .rc_condition_pando_projection

.rc_condition_pando_projection <- function(
    grn_result, membership, unit_meta, genes, comparison_support) {
  answer <- .rc_condition_pando_projection_strict_base(
    grn_result = grn_result,
    membership = membership,
    unit_meta = unit_meta,
    genes = genes,
    comparison_support = comparison_support
  )
  reliability <- matrix(
    NA_real_, nrow(answer$primary), ncol(answer$primary),
    dimnames = dimnames(answer$primary)
  )
  for (fit in grn_result$condition_grn_fits) {
    .rc_require_pando_condition_grn_fit(fit)
    coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
    condition_col <- as.character(fit$condition_col)[[1L]]
    celltype_col <- as.character(fit$cell_type_col)[[1L]]
    if (!all(c(condition_col, celltype_col) %in% colnames(unit_meta))) {
      stop("Metacell metadata lack fitted condition or cell-type columns.",
           call. = FALSE)
    }
    for (condition in fit$condition_levels) {
      one <- coefficient[
        as.character(coefficient$condition) == condition &
          coefficient$estimable %in% TRUE &
          is.finite(as.numeric(coefficient$padj)) &
          as.numeric(coefficient$padj) < 0.05,
        , drop = FALSE
      ]
      targets <- intersect(
        tolower(unique(as.character(one$target))), rownames(reliability)
      )
      units <- as.character(unit_meta$unit_id)[
        as.character(unit_meta[[condition_col]]) == condition &
          as.character(unit_meta[[celltype_col]]) == fit$cell_type
      ]
      if (length(targets) && length(units)) {
        reliability[targets, units] <- 1
      }
    }
  }
  answer$reliability <- reliability
  answer$origin <- "paired_cell_full_fit_fixed_dictionary_glm_padj_filtered"
  answer$primary_projection <- "padj_filtered_fixed_dictionary_condition_glm"
  answer$common_projection_role <-
    "compatibility_alias_of_primary_no_common_support_decomposition"
  answer$nonestimable_projection_policy <-
    "coefficient_NA_and_zero_realized_penalty_contribution"
  answer
}
