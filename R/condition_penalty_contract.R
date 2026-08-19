# Pando edge eligibility used by the canonical Stage-1 merge.
#
# Conditional GRNs use Pando Scheme E z=0.25 on one frozen exact-edge
# dictionary. Condition-wise BH and target R2 remain diagnostics only; neither
# changes conditional network membership or overwrites a continuous coefficient.
# Standard (single-condition) Pando keeps its existing BH/R2 filter below.

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
  if (!is.data.frame(fit_table) ||
      !all(required %in% colnames(fit_table)) ||
      !all(c("target", "condition") %in% colnames(coefficient))) {
    stop(
      "Condition-GRN diagnostics require target-level fit_status and full-data ",
      "rsq aligned to coefficient target and condition.", call. = FALSE
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
      "Condition-GRN coefficients cannot be aligned to target-level fit diagnostics.",
      call. = FALSE
    )
  }
  status <- trimws(as.character(fit_table$fit_status[index]))
  rsq <- suppressWarnings(as.numeric(fit_table$rsq[index]))
  if (anyNA(status) || any(!nzchar(status))) {
    stop("Condition-GRN fit_status values must be complete.", call. = FALSE)
  }
  data.frame(
    fit_status = status,
    target_rsq = rsq,
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
    "estimate", "estimable", "padj", "statistically_supported",
    "global_support", "local_support", "active", "significant",
    "penalty_effect", "contrast_identifiable", "shared_by_boundary",
    "fused_by_penalty", "penalty_family", "penalty_value",
    "solver_status", "kkt_residual", "iterations"
  )
  if (!is.data.frame(coefficient) ||
      !all(required %in% colnames(coefficient))) {
    stop(
      "Condition-GRN coefficients must retain the Scheme-E continuous-effect, ",
      "diagnostic-inference, identifiability and solver contract.",
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
  expected_statistical <- is.finite(estimate) & is.finite(padj) &
    padj < threshold
  expected_active <- is.finite(estimate)
  if (!identical(
      as.logical(coefficient$statistically_supported), expected_statistical
  )) {
    stop(
      "Pando statistically_supported flags do not match diagnostic condition-wise BH values.",
      call. = FALSE
    )
  }
  if (!identical(as.logical(coefficient$active), expected_active)) {
    stop(
      "Pando Scheme-E active flags must retain every finite coefficient on the ",
      "frozen common dictionary.", call. = FALSE
    )
  }
  if (!identical(
      as.logical(coefficient$significant), expected_statistical
  )) {
    stop(
      "Pando significant flags must remain diagnostic BH annotations and must ",
      "not define Scheme-E network membership.", call. = FALSE
    )
  }
  observed_effect <- suppressWarnings(as.numeric(coefficient$penalty_effect))
  expected_effect <- ifelse(expected_active, estimate, NA_real_)
  comparable <- is.finite(expected_effect) & is.finite(observed_effect)
  if (any(is.finite(expected_effect) != is.finite(observed_effect)) ||
      any(abs(expected_effect[comparable] - observed_effect[comparable]) > 1e-12)) {
    stop(
      "Pando penalty_effect must equal the continuous Scheme-E condition coefficient.",
      call. = FALSE
    )
  }
  family <- as.character(coefficient$penalty_family)
  value <- suppressWarnings(as.numeric(coefficient$penalty_value))
  solver <- as.character(coefficient$solver_status)
  kkt <- suppressWarnings(as.numeric(coefficient$kkt_residual))
  iteration <- suppressWarnings(as.integer(coefficient$iterations))
  if (anyNA(family) || any(family != "exact_edge_sparse_deviation") ||
      any(!is.finite(value)) || any(abs(value - 0.25) > 1e-15) ||
      anyNA(solver) || any(solver != "ok") ||
      any(!is.finite(kkt)) || any(kkt < 0) ||
      anyNA(iteration) || any(iteration < 0L)) {
    stop("Pando conditional coefficients are not a converged fixed z=0.25 Scheme-E fit.",
         call. = FALSE)
  }
  boundary <- coefficient$shared_by_boundary %in% TRUE
  identifiable <- coefficient$contrast_identifiable %in% TRUE
  if (any(boundary & identifiable) || any(!boundary & !identifiable)) {
    stop(
      "Scheme-E shared_by_boundary must be the complement of contrast_identifiable.",
      call. = FALSE
    )
  }
  invisible(expected_active)
}

.rc_condition_target_rsq <- function(coefficient) {
  if ("target_rsq" %in% colnames(coefficient)) {
    return(suppressWarnings(as.numeric(coefficient$target_rsq)))
  }
  if ("rsq" %in% colnames(coefficient)) {
    return(suppressWarnings(as.numeric(coefficient$rsq)))
  }
  stop(
    "Condition-GRN diagnostics require full-data target R2 (`rsq`).",
    call. = FALSE
  )
}

.rc_condition_penalty_gate <- function(
    coefficient, padj_threshold = NULL, target_rsq_threshold = NULL) {
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
  # Keep the threshold validated and exported as a diagnostic contract, but do
  # not use it to remove a Scheme-E common-dictionary edge.
  invisible(.rc_target_rsq_threshold(
    target_rsq_threshold %||% getOption(
      "RegCompassR.target_rsq_threshold", .RC_PANDO_TARGET_RSQ_THRESHOLD
    )
  ))
  active <- .rc_validate_pando_active_condition_edges(
    coefficient, padj_threshold = threshold
  )
  fit_status <- if ("fit_status" %in% colnames(coefficient)) {
    trimws(as.character(coefficient$fit_status))
  } else {
    rep("ok", nrow(coefficient))
  }
  effect <- suppressWarnings(as.numeric(coefficient$penalty_effect))
  active & !is.na(fit_status) & fit_status == "ok" & is.finite(effect)
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
  # Diagnostic only: this flag can be inspected, plotted or reported, but it is
  # not part of the Scheme-E edge/projection gate.
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
  fit$regcompass_penalty_filter <-
    "finite continuous Scheme-E coefficient on frozen dictionary & fit_status == 'ok'"
  fit$regcompass_fit_status_filter <- "fit_status == 'ok'"
  fit$regcompass_target_rsq_filter <- paste0(
    "diagnostic only: rsq >= ", format(rsq_threshold, trim = TRUE)
  )
  fit$regcompass_target_rsq_definition <-
    "scheme_e_z025_full_data_R2_diagnostic"
  fit$regcompass_rank_deficient_policy <- paste(
    "non-identifiable condition contrast is exact-shared and flagged;",
    "continuous condition coefficient remains on the fixed dictionary"
  )
  fit$regcompass_significance_role <-
    "BH_and_target_R2_are_diagnostics_only;no_edge_reselection"
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
         call. = FALSE
    )
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
