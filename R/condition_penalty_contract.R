# Pando edge eligibility used by the canonical Stage-1 merge.

# Retained only for backward-compatible provenance fields. They are not applied
# as post-fit thresholds; zero records that no additional correlation or
# coefficient-size gate is imposed after Pando candidate screening and BH.
.RC_PANDO_PENALTY_CORR_THRESHOLD <- 0
.RC_PANDO_PENALTY_ESTIMATE_THRESHOLD <- 0

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
  coefficient$estimable %in% TRUE &
    is.finite(padj) & padj < 0.05 &
    is.finite(estimate)
}

.rc_apply_condition_penalty_gate <- function(fit) {
  .rc_require_pando_condition_grn_fit(fit)
  coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
  gate <- .rc_condition_penalty_gate(coefficient)
  coefficient$significant <- gate
  coefficient$penalty_effect <- ifelse(
    gate, suppressWarnings(as.numeric(coefficient$estimate)), 0
  )
  fit$coefficients <- coefficient
  fit$regcompass_penalty_filter <- "estimable & BH padj < 0.05"
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
