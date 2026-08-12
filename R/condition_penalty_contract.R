# Pando edge eligibility used by the canonical Stage-1 merge.

.RC_PANDO_PENALTY_CORR_THRESHOLD <- 0
.RC_PANDO_PENALTY_ESTIMATE_THRESHOLD <- 0

.rc_condition_fit_status_for_coefficients <- function(fit, coefficient) {
  fit_table <- as.data.frame(fit$fit, stringsAsFactors = FALSE)
  required <- c("target", "condition", "fit_status")
  if (!is.data.frame(fit_table) ||
      !all(required %in% colnames(fit_table)) ||
      !all(c("target", "condition") %in% colnames(coefficient))) {
    stop(
      "Condition-GRN penalty filtering requires target-level fit_status ",
      "diagnostics aligned to coefficient target and condition.",
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
  if (anyNA(fit_key) || any(!nzchar(fit_key)) || anyDuplicated(fit_key)) {
    stop(
      "Condition-GRN target-level fit diagnostics are duplicated or incomplete.",
      call. = FALSE
    )
  }
  index <- match(coefficient_key, fit_key)
  if (anyNA(index)) {
    stop(
      "Condition-GRN coefficients cannot be aligned to target-level fit_status.",
      call. = FALSE
    )
  }
  status <- trimws(as.character(fit_table$fit_status[index]))
  if (anyNA(status) || any(!nzchar(status))) {
    stop("Condition-GRN fit_status values must be complete.", call. = FALSE)
  }
  status
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
    "estimate", "estimable", "padj", "statistically_supported",
    "global_support", "local_support", "active", "significant",
    "penalty_effect"
  )
  if (!is.data.frame(coefficient) ||
      !all(required %in% colnames(coefficient))) {
    stop(
      "Condition-GRN coefficients must retain Pando ridge statistics, ",
      "global/local candidate support, active flags, and penalty_effect.",
      call. = FALSE
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
  padj <- suppressWarnings(as.numeric(coefficient$padj))
  expected_statistical <- coefficient$estimable %in% TRUE &
    is.finite(estimate) & is.finite(padj) & padj < threshold
  expected_dictionary_support <-
    coefficient$global_support %in% TRUE |
    coefficient$local_support %in% TRUE
  expected_active <- expected_statistical & expected_dictionary_support
  if (!identical(
      as.logical(coefficient$statistically_supported), expected_statistical
  )) {
    stop(
      "Pando statistically_supported flags do not match condition-wise BH ridge evidence.",
      call. = FALSE
    )
  }
  if (!identical(as.logical(coefficient$active), expected_active) ||
      !identical(as.logical(coefficient$significant), expected_active)) {
    stop(
      "Pando active condition-edge flags must equal BH-supported ridge evidence ",
      "with pooled/global or condition-local Pando candidate support.",
      call. = FALSE
    )
  }
  expected_effect <- ifelse(expected_active, estimate, 0)
  observed_effect <- suppressWarnings(as.numeric(coefficient$penalty_effect))
  comparable <- is.finite(expected_effect) & is.finite(observed_effect)
  if (any(is.finite(expected_effect) != is.finite(observed_effect)) ||
      any(abs(expected_effect[comparable] - observed_effect[comparable]) > 1e-12)) {
    stop(
      "Pando penalty_effect must equal the active condition-specific ridge ",
      "coefficient or zero.", call. = FALSE
    )
  }
  invisible(expected_active)
}

.rc_condition_penalty_gate <- function(coefficient, padj_threshold = NULL) {
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
  active <- .rc_validate_pando_active_condition_edges(
    coefficient, padj_threshold = threshold
  )
  fit_status <- if ("fit_status" %in% colnames(coefficient)) {
    trimws(as.character(coefficient$fit_status))
  } else {
    rep("ok", nrow(coefficient))
  }
  active & !is.na(fit_status) & fit_status == "ok"
}

.rc_apply_condition_penalty_gate <- function(fit) {
  .rc_require_pando_condition_grn_fit(fit)
  threshold <- .rc_condition_padj_threshold(fit = fit)
  coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
  coefficient$fit_status <- .rc_condition_fit_status_for_coefficients(
    fit, coefficient
  )
  coefficient$padj_threshold <- threshold
  gate <- .rc_condition_penalty_gate(
    coefficient, padj_threshold = threshold
  )
  coefficient$penalty_eligible <- gate
  coefficient$active_in_condition <- gate
  # Preserve Pando's active/significant/penalty_effect contract. RegCompass adds
  # only the target-level fit_status == 'ok' requirement for downstream use.
  fit$coefficients <- coefficient
  fit$regcompass_penalty_filter <-
    "Pando active edge & target fit_status == 'ok'"
  fit$regcompass_fit_status_filter <- "fit_status == 'ok'"
  fit$regcompass_rank_deficient_policy <-
    "regularized_ok_fit_retained; non-estimable condition edge excluded"
  fit$regcompass_significance_role <-
    "consume_pando_active_condition_edge_without_reselection"
  fit$regcompass_padj_threshold <- threshold
  fit
}

.rc_filter_standard_pando_edges <- function(
    table, padj_threshold = .rc_standard_pando_padj_default) {
  required <- c("estimate", "padj")
  if (!is.data.frame(table) || !all(required %in% colnames(table))) {
    stop(
      "Standard Pando requires estimate and padj columns for RegCompass ",
      "penalty filtering.", call. = FALSE
    )
  }
  threshold <- suppressWarnings(as.numeric(padj_threshold))
  if (length(threshold) != 1L || !is.finite(threshold) ||
      threshold <= 0 || threshold >= 1) {
    stop("Standard Pando `padj_threshold` must be one value in (0, 1).",
         call. = FALSE)
  }
  estimate <- suppressWarnings(as.numeric(table$estimate))
  padj <- suppressWarnings(as.numeric(table$padj))
  keep <- is.finite(estimate) & is.finite(padj) & padj < threshold
  if ("estimable" %in% colnames(table)) {
    keep <- keep & table$estimable %in% TRUE
  }
  answer <- table[keep, , drop = FALSE]
  attr(answer, "edge_filter") <- list(
    estimable = if ("estimable" %in% colnames(table)) TRUE else NA,
    padj = paste0("< ", format(threshold, trim = TRUE))
  )
  answer
}
