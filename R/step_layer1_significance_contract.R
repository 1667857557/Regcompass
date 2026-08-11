# Significance-aware Layer-1 contract for condition multi-task ridge fits.
# Loaded after step_layer1_parallel_object_contract.R so the canonical
# condition projection keeps beta * mean(TF) * mean(ATAC) while edge inclusion
# and target reliability follow the configured BH threshold.

.rc_require_layer1_condition_grn_fit <- function(fit) {
  .rc_require_pando_condition_grn_fit(fit)
  if (is.null(fit$regcompass_penalty_filter)) return(invisible(TRUE))

  threshold <- .rc_condition_padj_threshold(fit = fit)
  expected_filter <- paste0(
    "estimable & finite estimate & fit_status == 'ok' & BH padj < ",
    format(threshold, trim = TRUE)
  )
  filter_value <- as.character(fit$regcompass_penalty_filter)
  if (length(filter_value) != 1L || is.na(filter_value) ||
      !identical(filter_value, expected_filter)) {
    stop(
      "RegCompass condition-GRN penalty gate metadata are inconsistent.",
      call. = FALSE
    )
  }

  coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
  required <- c(
    "significant", "penalty_effect", "estimate", "estimable", "padj"
  )
  if (!all(required %in% colnames(coefficient))) {
    stop(
      "RegCompass-gated condition fits require ridge coefficients, padj and ",
      "projection effects.", call. = FALSE
    )
  }
  expected_gate <- .rc_condition_penalty_gate(
    coefficient, padj_threshold = threshold
  )
  if (!identical(as.logical(coefficient$significant), expected_gate)) {
    stop("RegCompass significant flags do not match the final BH gate.",
         call. = FALSE)
  }
  estimate <- suppressWarnings(as.numeric(coefficient$estimate))
  expected_effect <- ifelse(expected_gate, estimate, 0)
  observed_effect <- suppressWarnings(as.numeric(coefficient$penalty_effect))
  comparable <- is.finite(expected_effect) & is.finite(observed_effect)
  if (any(is.finite(expected_effect) != is.finite(observed_effect)) ||
      any(abs(expected_effect[comparable] - observed_effect[comparable]) > 1e-12)) {
    stop(
      "RegCompass-gated penalty_effect does not match the final BH-significant ",
      "ridge gate.", call. = FALSE
    )
  }
  invisible(TRUE)
}

.rc_condition_target_reliability <- function(fit) {
  fit_table <- as.data.frame(fit$fit, stringsAsFactors = FALSE)
  coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
  required_fit <- c("target", "condition", "rsq", "fit_status")
  required_coefficient <- c(
    "target", "condition", "significant", "estimable", "estimate", "padj"
  )
  if (!all(required_fit %in% colnames(fit_table)) || !nrow(fit_table) ||
      !all(required_coefficient %in% colnames(coefficient))) {
    stop(
      "Condition-GRN reliability requires target-condition fit diagnostics ",
      "and BH-significant ridge-edge flags.", call. = FALSE
    )
  }

  fit_target <- toupper(trimws(as.character(fit_table$target)))
  fit_condition <- as.character(fit_table$condition)
  fit_key <- paste(fit_target, fit_condition, sep = "\001")
  if (anyNA(fit_key) || any(!nzchar(fit_target)) ||
      any(!nzchar(fit_condition)) || anyDuplicated(fit_key)) {
    stop(
      "Condition-GRN target-condition fit diagnostics must be complete and unique.",
      call. = FALSE
    )
  }

  coefficient_target <- toupper(trimws(as.character(coefficient$target)))
  coefficient_condition <- as.character(coefficient$condition)
  coefficient_key <- paste(
    coefficient_target, coefficient_condition, sep = "\001"
  )
  coefficient_index <- match(coefficient_key, fit_key)
  if (anyNA(coefficient_index)) {
    stop(
      "Condition-GRN coefficients cannot be aligned to target fit diagnostics.",
      call. = FALSE
    )
  }

  threshold <- .rc_condition_padj_threshold(fit = fit)
  final_gate <- .rc_condition_penalty_gate(
    transform(
      coefficient,
      fit_status = fit_table$fit_status[coefficient_index],
      padj_threshold = threshold
    ),
    padj_threshold = threshold
  )
  n_significant_edges <- tabulate(
    coefficient_index[final_gate], nbins = nrow(fit_table)
  )
  rsq <- suppressWarnings(as.numeric(fit_table$rsq))
  fit_status <- trimws(as.character(fit_table$fit_status))
  if (anyNA(fit_status) || any(!nzchar(fit_status))) {
    stop("Condition-GRN fit_status values must be complete.", call. = FALSE)
  }
  reliability <- rep(NA_real_, nrow(fit_table))
  eligible <- fit_status == "ok" &
    n_significant_edges > 0L & is.finite(rsq)
  reliability[eligible] <- sqrt(pmin(1, pmax(0, rsq[eligible])))

  data.frame(
    target = as.character(fit_table$target),
    condition = fit_condition,
    rsq = rsq,
    fit_status = fit_status,
    n_significant_edges = as.integer(n_significant_edges),
    n_projection_edges = as.integer(n_significant_edges),
    padj_threshold = threshold,
    reliability = reliability,
    stringsAsFactors = FALSE
  )
}

.rc_condition_pando_projection_bh_base <- .rc_condition_pando_projection

.rc_condition_pando_projection <- function(
    grn_result, membership, unit_meta, genes, rna_assay, atac_assay) {
  answer <- .rc_condition_pando_projection_bh_base(
    grn_result = grn_result,
    membership = membership,
    unit_meta = unit_meta,
    genes = genes,
    rna_assay = rna_assay,
    atac_assay = atac_assay
  )

  if (is.data.frame(answer$coverage) && nrow(answer$coverage)) {
    thresholds <- stats::setNames(vapply(
      grn_result$condition_grn_fits,
      function(fit) .rc_condition_padj_threshold(fit = fit),
      numeric(1)
    ), vapply(
      grn_result$condition_grn_fits,
      function(fit) as.character(fit$cell_type)[[1L]],
      character(1)
    ))
    answer$coverage$padj_threshold <- unname(
      thresholds[as.character(answer$coverage$cell_type)]
    )
    answer$coverage$n_significant_edges <-
      as.integer(answer$coverage$n_projection_edges)
    if ("padj_threshold_diagnostic" %in% colnames(answer$coverage)) {
      answer$coverage$padj_threshold_diagnostic <- NULL
    }
    if ("n_significant_diagnostic_edges" %in% colnames(answer$coverage)) {
      answer$coverage$n_significant_diagnostic_edges <- NULL
    }
    answer$coverage$projection_effect <-
      "BH_significant_multitask_ridge_penalty_effect"
  }
  answer$origin <-
    "paired_cell_significant_union_multitask_ridge_bh_filtered"
  answer$projection_name <-
    "bh_significant_multitask_ridge_condition_effect"
  answer$nonestimable_policy <-
    "nonestimable_or_nonsignificant_condition_edge_has_zero_projection_contribution"
  answer
}
