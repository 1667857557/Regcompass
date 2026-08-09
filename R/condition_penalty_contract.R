# Pando edge eligibility used by the canonical Stage-1 merge.

# Retained only for backward-compatible provenance fields. They are not applied
# as post-fit thresholds; zero records that no additional correlation or
# coefficient-size gate is imposed after Pando candidate screening and BH.
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

.rc_condition_penalty_gate <- function(coefficient) {
  required <- c("estimate", "padj", "estimable")
  if (!is.data.frame(coefficient) ||
      !all(required %in% colnames(coefficient))) {
    stop(
      "Condition-GRN coefficients must contain estimate, padj, and estimable ",
      "columns before RegCompass penalty filtering.",
      call. = FALSE
    )
  }
  estimate <- suppressWarnings(as.numeric(coefficient$estimate))
  padj <- suppressWarnings(as.numeric(coefficient$padj))

  # Canonical condition-GRN fits attach the target-level status before this
  # gate is applied. The fallback preserves compatibility for legacy tables
  # that predate fit_status metadata; it is not used by the current condition
  # merge/projection path.
  fit_status <- if ("fit_status" %in% colnames(coefficient)) {
    trimws(as.character(coefficient$fit_status))
  } else {
    rep("ok", nrow(coefficient))
  }

  coefficient$estimable %in% TRUE &
    !is.na(fit_status) & fit_status == "ok" &
    is.finite(padj) & padj < 0.05 &
    is.finite(estimate)
}

.rc_apply_condition_penalty_gate <- function(fit) {
  .rc_require_pando_condition_grn_fit(fit)
  coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
  coefficient$fit_status <- .rc_condition_fit_status_for_coefficients(
    fit, coefficient
  )
  gate <- .rc_condition_penalty_gate(coefficient)
  coefficient$significant <- gate
  coefficient$penalty_effect <- ifelse(
    gate, suppressWarnings(as.numeric(coefficient$estimate)), 0
  )
  fit$coefficients <- coefficient
  # Keep the historical BH-filter field for compatibility and record the new
  # target-level eligibility rule separately so older saved objects remain
  # readable while current fits are explicit about rank-deficient exclusion.
  fit$regcompass_penalty_filter <- "estimable & BH padj < 0.05"
  fit$regcompass_fit_status_filter <- "fit_status == 'ok'"
  fit$regcompass_rank_deficient_policy <- "exclude_from_penalty"
  fit
}

.rc_filter_standard_pando_edges <- function(table) {
  required <- c("estimate", "padj")
  if (!is.data.frame(table) || !all(required %in% colnames(table))) {
    stop(
      "Standard Pando requires estimate and padj columns for RegCompass ",
      "penalty filtering.",
      call. = FALSE
    )
  }
  estimate <- suppressWarnings(as.numeric(table$estimate))
  padj <- suppressWarnings(as.numeric(table$padj))
  keep <- is.finite(estimate) &
    is.finite(padj) & padj < .rc_standard_pando_padj_fixed
  if ("estimable" %in% colnames(table)) {
    keep <- keep & table$estimable %in% TRUE
  }
  answer <- table[keep, , drop = FALSE]
  attr(answer, "edge_filter") <- list(
    estimable = if ("estimable" %in% colnames(table)) TRUE else NA,
    padj = "< 0.05"
  )
  answer
}
