# Tighten active-edge semantics after the core multitask fitter is loaded.

.rc_fit_multitask_target_core <- .rc_fit_multitask_target

.rc_fit_multitask_target <- function(
    edges, target, rna, atac, meta, condition_col, args) {
  answer <- .rc_fit_multitask_target_core(
    edges = edges,
    target = target,
    rna = rna,
    atac = atac,
    meta = meta,
    condition_col = condition_col,
    args = args
  )
  condition <- answer$condition
  if (!is.data.frame(condition) || !nrow(condition)) return(answer)

  effect_threshold <- max(
    as.numeric(args$min_abs_effect),
    as.numeric(args$zero_tolerance)
  )
  bootstrap_fraction <- suppressWarnings(
    as.numeric(condition$n_bootstrap_success) /
      as.numeric(condition$n_bootstrap_requested)
  )
  bootstrap_adequate <- is.finite(bootstrap_fraction) &
    bootstrap_fraction >= args$min_bootstrap_success_fraction
  condition$bootstrap_success_fraction <- bootstrap_fraction
  condition$bootstrap_completion_adequate <- bootstrap_adequate
  condition$active_edge <-
    bootstrap_adequate &
    is.finite(condition$selection_frequency) &
    condition$selection_frequency >= args$min_selection_frequency &
    is.finite(condition$sign_stability) &
    condition$sign_stability >= args$min_sign_stability &
    is.finite(condition$effective_estimate) &
    abs(condition$effective_estimate) > effect_threshold &
    is.finite(condition$cv_rsq) &
    condition$cv_rsq >= args$min_cv_rsq
  condition$stable_estimate[
    !is.finite(condition$stable_estimate)
  ] <- 0
  condition$estimate <- condition$stable_estimate

  sign_by_edge <- split(condition, condition$edge_id)
  flip <- vapply(sign_by_edge, function(one) {
    one <- one[one$active_edge %in% TRUE, , drop = FALSE]
    length(unique(sign(one$effective_estimate))) > 1L
  }, logical(1))
  condition$sign_flip_flag <- unname(flip[condition$edge_id])
  answer$condition <- condition

  global <- answer$global
  if (is.data.frame(global) && nrow(global)) {
    global$bootstrap_success_fraction <- bootstrap_fraction[[1L]]
    global$bootstrap_completion_adequate <- bootstrap_adequate[[1L]]
    answer$global <- global
  }

  diagnostics <- answer$diagnostics
  if (is.data.frame(diagnostics) && nrow(diagnostics)) {
    diagnostics$n_active_condition_edges <- sum(
      condition$active_edge %in% TRUE
    )
    diagnostics$n_active_conditions <- length(unique(
      condition$condition[condition$active_edge %in% TRUE]
    ))
    diagnostics$active_effect_threshold <- effect_threshold
    diagnostics$bootstrap_success_fraction <- bootstrap_fraction[[1L]]
    diagnostics$min_bootstrap_success_fraction <-
      args$min_bootstrap_success_fraction
    diagnostics$bootstrap_completion_adequate <- bootstrap_adequate[[1L]]
    answer$diagnostics <- diagnostics
  }
  answer
}
