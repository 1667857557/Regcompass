# Pando edge eligibility used by the canonical Stage-1 merge.

.RC_PANDO_PENALTY_CORR_THRESHOLD <- 0
.RC_PANDO_PENALTY_ESTIMATE_THRESHOLD <- 0
.RC_PANDO_TARGET_RSQ_THRESHOLD_DEFAULT <- 0.05

.rc_target_rsq_threshold <- function(value = .RC_PANDO_TARGET_RSQ_THRESHOLD_DEFAULT) {
  value <- suppressWarnings(as.numeric(value))
  if (length(value) != 1L || !is.finite(value) || value < 0 || value > 1) {
    stop("`target_rsq_threshold` must be one finite value in [0, 1].",
         call. = FALSE)
  }
  value[[1L]]
}

.rc_condition_fit_diagnostics_for_coefficients <- function(
    fit, coefficient, target_rsq_threshold = .RC_PANDO_TARGET_RSQ_THRESHOLD_DEFAULT) {
  threshold <- .rc_target_rsq_threshold(target_rsq_threshold)
  fit_table <- as.data.frame(fit$fit, stringsAsFactors = FALSE)
  required <- c(
    "target", "condition", "fit_status", "rsq", "rsq_in_sample", "rsq_oof"
  )
  if (!is.data.frame(fit_table) || !all(required %in% colnames(fit_table)) ||
      !all(c("target", "condition") %in% colnames(coefficient))) {
    stop(
      "Condition-GRN penalty filtering requires target-level full-data and ",
      "OOF fit diagnostics aligned to coefficient target and condition.",
      call. = FALSE
    )
  }
  rsq <- suppressWarnings(as.numeric(fit_table$rsq))
  rsq_in_sample <- suppressWarnings(as.numeric(fit_table$rsq_in_sample))
  comparable <- is.finite(rsq) & is.finite(rsq_in_sample)
  if (any(is.finite(rsq) != is.finite(rsq_in_sample)) ||
      any(abs(rsq[comparable] - rsq_in_sample[comparable]) > 1e-12)) {
    stop(
      "Pando `rsq` must equal the selected-lambda final full-data R-squared; ",
      "install the matching Pando_regcompass ridge-diagnostics revision.",
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
      "Condition-GRN coefficients cannot be aligned to target-level fit diagnostics.",
      call. = FALSE
    )
  }
  status <- trimws(as.character(fit_table$fit_status[index]))
  target_rsq <- rsq[index]
  target_rsq_oof <- suppressWarnings(as.numeric(fit_table$rsq_oof[index]))
  if (anyNA(status) || any(!nzchar(status))) {
    stop("Condition-GRN fit_status values must be complete.", call. = FALSE)
  }
  # "evaluated" records that a final target model was successfully fit. A
  # non-finite R2 is therefore evaluated-but-unsupported rather than unavailable.
  evaluated <- status == "ok"
  supported <- evaluated & is.finite(target_rsq) & target_rsq >= threshold
  data.frame(
    fit_status = status,
    target_rsq = target_rsq,
    target_rsq_oof = target_rsq_oof,
    target_model_evaluated = evaluated,
    target_model_supported = supported,
    stringsAsFactors = FALSE
  )
}

.rc_condition_fit_status_for_coefficients <- function(fit, coefficient) {
  .rc_condition_fit_diagnostics_for_coefficients(
    fit, coefficient,
    target_rsq_threshold = .RC_PANDO_TARGET_RSQ_THRESHOLD_DEFAULT
  )$fit_status
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
      "candidate-support provenance, active flags, and penalty_effect.",
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
  expected_active <- expected_statistical
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
      "Pando active condition-edge flags must equal the condition's own ",
      "estimable BH-supported ridge evidence; global/local correlation support ",
      "is candidate provenance only.",
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

# Edge-level Pando/fit-status gate only. Target-model R2 is intentionally kept
# separate so RegCompass does not redefine Pando's active/significant statistic.
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
    rep(NA_character_, nrow(coefficient))
  }
  active & !is.na(fit_status) & fit_status == "ok"
}

.rc_apply_condition_penalty_gate <- function(
    fit, target_rsq_threshold = .RC_PANDO_TARGET_RSQ_THRESHOLD_DEFAULT) {
  .rc_require_pando_condition_grn_fit(fit)
  threshold <- .rc_condition_padj_threshold(fit = fit)
  target_threshold <- .rc_target_rsq_threshold(target_rsq_threshold)
  coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
  diagnostics <- .rc_condition_fit_diagnostics_for_coefficients(
    fit, coefficient, target_rsq_threshold = target_threshold
  )
  coefficient$fit_status <- diagnostics$fit_status
  coefficient$target_rsq <- diagnostics$target_rsq
  coefficient$target_rsq_oof <- diagnostics$target_rsq_oof
  coefficient$target_model_evaluated <- diagnostics$target_model_evaluated
  coefficient$target_model_supported <- diagnostics$target_model_supported
  coefficient$padj_threshold <- threshold
  coefficient$target_rsq_threshold <- target_threshold
  edge_gate <- .rc_condition_penalty_gate(
    coefficient, padj_threshold = threshold
  )
  gate <- edge_gate & coefficient$target_model_supported %in% TRUE
  coefficient$penalty_eligible <- gate
  coefficient$active_in_condition <- gate
  fit$coefficients <- coefficient
  fit$regcompass_penalty_filter <- paste0(
    "Pando BH-active edge & fit_status == 'ok' & final full-data R2 >= ",
    format(target_threshold, trim = TRUE)
  )
  fit$regcompass_fit_status_filter <- "fit_status == 'ok'"
  fit$regcompass_target_rsq_threshold <- target_threshold
  fit$regcompass_target_rsq_role <-
    "selected_lambda_final_full_data_target_model_quality_gate"
  fit$regcompass_oof_rsq_role <- "diagnostic_only"
  fit$regcompass_rank_deficient_policy <-
    "regularized_ok_fit_retained; non-estimable condition edge excluded"
  fit$regcompass_significance_role <-
    "consume_pando_condition_bh_active_edge_without_reselection"
  fit$regcompass_padj_threshold <- threshold
  fit
}

.rc_annotate_standard_pando_edges <- function(
    table, target_rsq_threshold = .RC_PANDO_TARGET_RSQ_THRESHOLD_DEFAULT) {
  target_threshold <- .rc_target_rsq_threshold(target_rsq_threshold)
  required <- c("estimate", "padj", "rsq", "rsq_in_sample")
  if (!is.data.frame(table) || !all(required %in% colnames(table))) {
    stop(
      "Standard Pando penalty filtering requires estimate, padj, rsq, and ",
      "rsq_in_sample columns.", call. = FALSE
    )
  }
  rsq <- suppressWarnings(as.numeric(table$rsq))
  rsq_in_sample <- suppressWarnings(as.numeric(table$rsq_in_sample))
  comparable <- is.finite(rsq) & is.finite(rsq_in_sample)
  if (any(is.finite(rsq) != is.finite(rsq_in_sample)) ||
      any(abs(rsq[comparable] - rsq_in_sample[comparable]) > 1e-12)) {
    stop(
      "Standard Pando `rsq` must equal selected-lambda final full-data R-squared.",
      call. = FALSE
    )
  }
  status <- if ("fit_status" %in% colnames(table)) {
    trimws(as.character(table$fit_status))
  } else {
    rep("ok", nrow(table))
  }
  table$target_rsq <- rsq
  table$target_rsq_oof <- if ("rsq_oof" %in% colnames(table)) {
    suppressWarnings(as.numeric(table$rsq_oof))
  } else {
    NA_real_
  }
  table$target_model_evaluated <- !is.na(status) & status == "ok"
  table$target_model_supported <- table$target_model_evaluated &
    is.finite(rsq) & rsq >= target_threshold
  table$target_rsq_threshold <- target_threshold
  table
}

.rc_filter_standard_pando_edges <- function(
    table, padj_threshold = .rc_standard_pando_padj_default,
    target_rsq_threshold = .RC_PANDO_TARGET_RSQ_THRESHOLD_DEFAULT) {
  threshold <- suppressWarnings(as.numeric(padj_threshold))
  if (length(threshold) != 1L || !is.finite(threshold) ||
      threshold <= 0 || threshold >= 1) {
    stop("Standard Pando `padj_threshold` must be one value in (0, 1).",
         call. = FALSE)
  }
  annotated <- .rc_annotate_standard_pando_edges(
    table, target_rsq_threshold = target_rsq_threshold
  )
  estimate <- suppressWarnings(as.numeric(annotated$estimate))
  padj <- suppressWarnings(as.numeric(annotated$padj))
  keep <- is.finite(estimate) & is.finite(padj) & padj < threshold &
    annotated$target_model_supported %in% TRUE
  if ("estimable" %in% colnames(annotated)) {
    keep <- keep & annotated$estimable %in% TRUE
  }
  annotated$penalty_eligible <- keep
  annotated$active_in_condition <- keep
  answer <- annotated[keep, , drop = FALSE]
  attr(answer, "edge_filter") <- list(
    estimable = if ("estimable" %in% colnames(annotated)) TRUE else NA,
    padj = paste0("< ", format(threshold, trim = TRUE)),
    target_rsq = paste0(
      ">= ", format(.rc_target_rsq_threshold(target_rsq_threshold), trim = TRUE)
    )
  )
  answer
}
