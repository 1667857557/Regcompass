# Binary target weight for condition-specific Pando regulatory effects.
#
# A target-condition pair contributes to the downstream regulatory modifier when
# it has at least one penalty-eligible active edge and the target fit completed
# successfully. Ridge R2/OOF R2 remain diagnostics only and never attenuate q.

.rc_condition_target_reliability <- function(fit) {
  fit_table <- as.data.frame(fit$fit, stringsAsFactors = FALSE)
  coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
  required_fit <- c("target", "condition", "rsq", "fit_status")
  required_coefficient <- c(
    "target", "condition", "active", "estimable", "estimate", "padj",
    "global_support", "local_support"
  )
  if (!all(required_fit %in% colnames(fit_table)) || !nrow(fit_table) ||
      !all(required_coefficient %in% colnames(coefficient))) {
    stop(
      "Condition-GRN target weight requires target-condition fit diagnostics ",
      "and Pando active ridge-edge flags.", call. = FALSE
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
  if (anyNA(coefficient_key) || any(!nzchar(coefficient_target)) ||
      any(!nzchar(coefficient_condition))) {
    stop("Condition-GRN coefficients contain incomplete target-condition labels.",
         call. = FALSE)
  }
  coefficient_index <- match(coefficient_key, fit_key)
  if (anyNA(coefficient_index)) {
    stop(
      "Condition-GRN coefficients cannot be aligned to target fit diagnostics.",
      call. = FALSE
    )
  }

  threshold <- .rc_condition_padj_threshold(fit = fit)
  coefficient$fit_status <- fit_table$fit_status[coefficient_index]
  coefficient$padj_threshold <- threshold
  final_gate <- .rc_condition_penalty_gate(
    coefficient, padj_threshold = threshold
  )
  n_active_edges <- tabulate(
    coefficient_index[final_gate], nbins = nrow(fit_table)
  )
  rsq <- suppressWarnings(as.numeric(fit_table$rsq))
  fit_status <- trimws(as.character(fit_table$fit_status))
  if (anyNA(fit_status) || any(!nzchar(fit_status))) {
    stop("Condition-GRN fit_status values must be complete.", call. = FALSE)
  }

  # Canonical target weight: q_{g,c}=1 for every valid active target.
  # OOF/in-sample R2 are retained below only for diagnostics and are not a gate.
  reliability <- rep(NA_real_, nrow(fit_table))
  eligible <- fit_status == "ok" & n_active_edges > 0L
  reliability[eligible] <- 1

  data.frame(
    target = as.character(fit_table$target),
    condition = fit_condition,
    rsq = rsq,
    fit_status = fit_status,
    n_significant_edges = as.integer(n_active_edges),
    n_projection_edges = as.integer(n_active_edges),
    padj_threshold = threshold,
    reliability = reliability,
    reliability_definition =
      "binary_active_target_q1; rsq_diagnostic_only",
    stringsAsFactors = FALSE
  )
}

# Keep the established projection implementation, but make its emitted coverage
# metadata agree with the binary q contract. Dynamic symbol lookup inside the
# core projection resolves .rc_condition_target_reliability() to the definition
# above, so the numeric downstream modifier also uses q=1.
.rc_condition_pando_projection_core <- .rc_condition_pando_projection
.rc_condition_pando_projection <- function(...) {
  out <- .rc_condition_pando_projection_core(...)
  if (is.data.frame(out$coverage) &&
      "reliability_definition" %in% colnames(out$coverage)) {
    out$coverage$reliability_definition <-
      "binary_active_target_q1; rsq_diagnostic_only"
  }
  out
}
