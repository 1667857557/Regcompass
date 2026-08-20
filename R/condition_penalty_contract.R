# Pando edge eligibility used by the canonical Stage-1 merge.
#
# Conditional GRNs use Pando E-star/JSE z=0.25 on one frozen exact-edge
# dictionary. Joint-refit Wald p-values are adjusted within condition x target.
# An exact edge enters RegCompass when any condition has padj < threshold and
# every condition has a valid production coefficient. Target R2 stays diagnostic.

.RC_PANDO_PENALTY_CORR_THRESHOLD <- 0
.RC_PANDO_PENALTY_ESTIMATE_THRESHOLD <- 0
.RC_PANDO_TARGET_RSQ_THRESHOLD <- 0.05

.rc_target_rsq_threshold <- function(value = getOption(
    "RegCompassR.target_rsq_threshold", .RC_PANDO_TARGET_RSQ_THRESHOLD)) {
  value <- suppressWarnings(as.numeric(value))
  if (length(value) != 1L || !is.finite(value) || value < 0 || value >= 1) {
    stop("RegCompass target R2 threshold must be one value in [0, 1).",
         call. = FALSE)
  }
  value[[1L]]
}

.rc_condition_fit_diagnostics_for_coefficients <- function(fit, coefficient) {
  fit_table <- as.data.frame(fit$fit, stringsAsFactors = FALSE)
  required <- c("target", "condition", "fit_status", "rsq")
  if (!all(required %in% colnames(fit_table)) ||
      !all(c("target", "condition") %in% colnames(coefficient))) {
    stop(
      "Condition-GRN diagnostics require target-level fit_status and rsq.",
      call. = FALSE
    )
  }
  fit_key <- paste(
    toupper(trimws(as.character(fit_table$target))),
    as.character(fit_table$condition), sep = "\001"
  )
  coefficient_key <- paste(
    toupper(trimws(as.character(coefficient$target))),
    as.character(coefficient$condition), sep = "\001"
  )
  if (anyNA(fit_key) || anyDuplicated(fit_key)) {
    stop("Condition-GRN target diagnostics are duplicated or incomplete.",
         call. = FALSE)
  }
  index <- match(coefficient_key, fit_key)
  if (anyNA(index)) {
    stop("Condition-GRN coefficients cannot be aligned to diagnostics.",
         call. = FALSE)
  }
  data.frame(
    fit_status = as.character(fit_table$fit_status[index]),
    target_rsq = suppressWarnings(as.numeric(fit_table$rsq[index])),
    stringsAsFactors = FALSE
  )
}

.rc_condition_padj_threshold <- function(fit = NULL, coefficient = NULL) {
  value <- if (!is.null(fit) && !is.null(fit$padj_threshold)) {
    fit$padj_threshold
  } else if (is.data.frame(coefficient) &&
             "padj_threshold" %in% colnames(coefficient)) {
    unique(suppressWarnings(as.numeric(coefficient$padj_threshold)))
  } else {
    0.05
  }
  value <- suppressWarnings(as.numeric(value))
  value <- value[is.finite(value)]
  if (length(value) != 1L || value <= 0 || value >= 1) {
    stop("Condition-GRN padj_threshold must be one value in (0, 1).",
         call. = FALSE)
  }
  value[[1L]]
}

.rc_validate_pando_active_condition_edges <- function(
    coefficient, padj_threshold = NULL) {
  required <- c(
    "edge_id", "target", "condition", "estimate", "penalty_effect",
    "pval", "padj", "inference_estimable", "condition_significant",
    "statistically_supported", "significant", "pando_estimation_active",
    "active", "edge_union_supported", "supporting_conditions",
    "n_supporting_conditions", "all_conditions_fit_valid",
    "active_in_regcompass", "fit_status", "penalty_family",
    "penalty_value", "solver_status", "kkt_residual", "iterations"
  )
  if (!is.data.frame(coefficient) ||
      !all(required %in% colnames(coefficient))) {
    stop(
      "Condition-GRN coefficients must retain the E-star/JSE inference and ",
      "exact-edge union handoff contract.", call. = FALSE
    )
  }
  threshold <- if (is.null(padj_threshold)) {
    .rc_condition_padj_threshold(coefficient = coefficient)
  } else {
    value <- suppressWarnings(as.numeric(padj_threshold))
    if (length(value) != 1L || !is.finite(value) ||
        value <= 0 || value >= 1) {
      stop("Condition-GRN padj_threshold must be in (0, 1).",
           call. = FALSE)
    }
    value
  }
  estimate <- suppressWarnings(as.numeric(coefficient$estimate))
  effect <- suppressWarnings(as.numeric(coefficient$penalty_effect))
  padj <- suppressWarnings(as.numeric(coefficient$padj))
  expected_condition_significant <-
    coefficient$inference_estimable %in% TRUE &
    is.finite(padj) & padj < threshold
  if (!identical(
      as.logical(coefficient$condition_significant),
      expected_condition_significant
  ) || !identical(
      as.logical(coefficient$statistically_supported),
      expected_condition_significant
  ) || !identical(
      as.logical(coefficient$significant),
      expected_condition_significant
  )) {
    stop("Condition significance does not match condition-target BH.",
         call. = FALSE)
  }
  expected_active <- is.finite(estimate)
  if (!identical(
      as.logical(coefficient$pando_estimation_active), expected_active
  ) || !identical(as.logical(coefficient$active), expected_active)) {
    stop("Pando estimation-active flags must track finite E-star coefficients.",
         call. = FALSE)
  }
  comparable <- is.finite(estimate) & is.finite(effect)
  if (any(is.finite(estimate) != is.finite(effect)) ||
      any(abs(estimate[comparable] - effect[comparable]) > 1e-12)) {
    stop("penalty_effect must equal the production E-star coefficient.",
         call. = FALSE)
  }

  edge_ids <- unique(as.character(coefficient$edge_id))
  for (edge_id in edge_ids) {
    index <- which(as.character(coefficient$edge_id) == edge_id)
    valid <- all(
      as.character(coefficient$fit_status[index]) == "ok" &
      is.finite(effect[index])
    )
    supporting <- as.character(coefficient$condition[index][
      expected_condition_significant[index]
    ])
    expected_union <- valid && length(supporting) > 0L
    if (any(as.logical(coefficient$all_conditions_fit_valid[index]) != valid) ||
        any(as.logical(coefficient$edge_union_supported[index]) != expected_union) ||
        any(as.logical(coefficient$active_in_regcompass[index]) != expected_union) ||
        any(as.integer(coefficient$n_supporting_conditions[index]) !=
            length(supporting))) {
      stop("Exact-edge RegCompass union flags are inconsistent.",
           call. = FALSE)
    }
    stored <- unique(as.character(coefficient$supporting_conditions[index]))
    expected_text <- paste(
      unique(as.character(coefficient$condition[index][
        expected_condition_significant[index]
      ])), collapse = ";"
    )
    if (length(stored) != 1L || !identical(stored, expected_text)) {
      stop("Exact-edge supporting-condition audit is inconsistent.",
           call. = FALSE)
    }
  }

  family <- as.character(coefficient$penalty_family)
  value <- suppressWarnings(as.numeric(coefficient$penalty_value))
  solver <- as.character(coefficient$solver_status)
  kkt <- suppressWarnings(as.numeric(coefficient$kkt_residual))
  iteration <- suppressWarnings(as.integer(coefficient$iterations))
  if (anyNA(family) ||
      any(family != .RC_PANDO_CONDITION_PENALTY_FAMILY) ||
      any(!is.finite(value)) ||
      any(abs(value - .RC_PANDO_CONDITION_SCHEME_E_Z) > 1e-15) ||
      anyNA(solver) || any(solver != "ok") ||
      any(!is.finite(kkt)) || any(kkt < 0) ||
      anyNA(iteration) || any(iteration < 0L)) {
    stop("Conditional coefficients are not converged E-star z=0.25 results.",
         call. = FALSE)
  }
  invisible(as.logical(coefficient$active_in_regcompass))
}

.rc_condition_target_rsq <- function(coefficient) {
  if ("target_rsq" %in% colnames(coefficient)) {
    return(suppressWarnings(as.numeric(coefficient$target_rsq)))
  }
  if ("rsq" %in% colnames(coefficient)) {
    return(suppressWarnings(as.numeric(coefficient$rsq)))
  }
  stop("Condition-GRN diagnostics require full-data target R2 (`rsq`).",
       call. = FALSE)
}

.rc_condition_penalty_gate <- function(
    coefficient, padj_threshold = NULL, target_rsq_threshold = NULL) {
  threshold <- if (is.null(padj_threshold)) {
    .rc_condition_padj_threshold(coefficient = coefficient)
  } else {
    .rc_condition_padj_threshold(
      coefficient = transform(coefficient, padj_threshold = padj_threshold)
    )
  }
  invisible(.rc_target_rsq_threshold(
    target_rsq_threshold %||% getOption(
      "RegCompassR.target_rsq_threshold", .RC_PANDO_TARGET_RSQ_THRESHOLD
    )
  ))
  gate <- .rc_validate_pando_active_condition_edges(
    coefficient, padj_threshold = threshold
  )
  effect <- suppressWarnings(as.numeric(coefficient$penalty_effect))
  gate & as.character(coefficient$fit_status) == "ok" & is.finite(effect)
}

.rc_apply_condition_penalty_gate <- function(
    fit, target_rsq_threshold = NULL) {
  .rc_require_pando_condition_grn_fit(fit)
  threshold <- .rc_condition_padj_threshold(fit = fit)
  rsq_threshold <- .rc_target_rsq_threshold(
    target_rsq_threshold %||% getOption(
      "RegCompassR.target_rsq_threshold", .RC_PANDO_TARGET_RSQ_THRESHOLD
    )
  )
  coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
  diagnostics <- .rc_condition_fit_diagnostics_for_coefficients(
    fit, coefficient
  )
  coefficient$fit_status <- diagnostics$fit_status
  coefficient$target_rsq <- diagnostics$target_rsq
  coefficient$rsq <- diagnostics$target_rsq
  coefficient$padj_threshold <- threshold
  coefficient$target_rsq_threshold <- rsq_threshold
  coefficient$target_model_supported <-
    coefficient$fit_status == "ok" &
    is.finite(coefficient$target_rsq) &
    coefficient$target_rsq >= rsq_threshold
  gate <- .rc_condition_penalty_gate(
    coefficient,
    padj_threshold = threshold,
    target_rsq_threshold = rsq_threshold
  )
  coefficient$penalty_eligible <- gate
  coefficient$active_in_condition <- gate
  fit$coefficients <- coefficient
  fit$regcompass_penalty_filter <- paste(
    "exact edge: all conditions valid and at least one condition has",
    paste0("condition-target BH padj < ", format(threshold, trim = TRUE))
  )
  fit$regcompass_fit_status_filter <- "fit_status == 'ok' in every condition"
  fit$regcompass_target_rsq_filter <- paste0(
    "diagnostic only: rsq >= ", format(rsq_threshold, trim = TRUE)
  )
  fit$regcompass_target_rsq_definition <-
    "scheme_e_z025_full_data_R2_diagnostic"
  fit$regcompass_significance_role <-
    "condition_target_BH_defines_any_condition_exact_edge_union;R2_diagnostic_only"
  fit$regcompass_padj_threshold <- threshold
  fit$regcompass_target_rsq_threshold <- rsq_threshold
  fit
}

.rc_filter_standard_pando_edges <- function(
    table, padj_threshold = .rc_standard_pando_padj_default,
    target_rsq_threshold = NULL) {
  required <- c("estimate", "padj", "rsq")
  if (!is.data.frame(table) || !all(required %in% colnames(table))) {
    stop(
      "Standard Pando requires estimate, padj and full-data rsq columns for ",
      "RegCompass penalty filtering.", call. = FALSE
    )
  }
  threshold <- suppressWarnings(as.numeric(padj_threshold))
  if (length(threshold) != 1L || !is.finite(threshold) ||
      threshold <= 0 || threshold >= 1) {
    stop("Standard Pando `padj_threshold` must be one value in (0, 1).",
         call. = FALSE)
  }
  rsq_threshold <- .rc_target_rsq_threshold(
    target_rsq_threshold %||% getOption(
      "RegCompassR.target_rsq_threshold", .RC_PANDO_TARGET_RSQ_THRESHOLD
    )
  )
  estimate <- suppressWarnings(as.numeric(table$estimate))
  padj <- suppressWarnings(as.numeric(table$padj))
  rsq <- suppressWarnings(as.numeric(table$rsq))
  keep <- is.finite(estimate) & is.finite(padj) & padj < threshold &
    is.finite(rsq) & rsq >= rsq_threshold
  if ("estimable" %in% colnames(table)) {
    keep <- keep & table$estimable %in% TRUE
  }
  answer <- table[keep, , drop = FALSE]
  attr(answer, "edge_filter") <- list(
    estimable = if ("estimable" %in% colnames(table)) TRUE else NA,
    padj = paste0("< ", format(threshold, trim = TRUE)),
    rsq = paste0(">= ", format(rsq_threshold, trim = TRUE)),
    rsq_definition = "selected_lambda_full_data_R2"
  )
  answer
}
