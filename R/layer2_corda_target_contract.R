# Apply the fixed original CORDA2 flux threshold to post-build scoring targets.

.rc_corda2_apply_target_flux <- function(
    model, flux_threshold = 1e-7, strict, cell_type) {
  diagnostics <- as.data.frame(model$closure_diagnostics)
  required <- c(
    "reaction_id", "target_direction", "feasible", "vmax",
    "final_feasible", "final_vmax"
  )
  if (!all(required %in% colnames(diagnostics))) {
    stop("CORDA2 closure diagnostics are incomplete.", call. = FALSE)
  }
  threshold <- as.numeric(flux_threshold)
  if (length(threshold) != 1L || !is.finite(threshold) || threshold <= 0) {
    stop("CORDA2 flux threshold must be positive.", call. = FALSE)
  }
  parent_target_feasible <- diagnostics$feasible %in% TRUE &
    is.finite(diagnostics$vmax) & diagnostics$vmax >= threshold
  final_target_feasible <- diagnostics$final_feasible %in% TRUE &
    is.finite(diagnostics$final_vmax) & diagnostics$final_vmax >= threshold
  failed <- parent_target_feasible & !final_target_feasible
  diagnostics$parent_corda2_feasible <- parent_target_feasible
  diagnostics$final_corda2_feasible <- final_target_feasible
  diagnostics$corda2_flux_threshold <- threshold
  diagnostics$completion_status <- ifelse(
    !parent_target_feasible,
    "parent_blocked",
    ifelse(final_target_feasible, "corda2_retained", "corda2_removed")
  )
  targets <- diagnostics[
    final_target_feasible,
    c("reaction_id", "target_direction"),
    drop = FALSE
  ]
  rownames(targets) <- NULL
  model$closure_diagnostics <- diagnostics
  model$target_directions <- targets
  model$target_status <- if (any(failed)) {
    "core_direction_removed_by_corda2"
  } else if (!nrow(targets)) {
    "no_feasible_core_direction"
  } else {
    "ok"
  }
  model$build_params$strict_requested <- strict
  model$build_params$strict_used_for_reconstruction <- FALSE
  model$build_params$n_parent_corda2_feasible_core_directions <-
    sum(parent_target_feasible)
  model$build_params$n_final_corda2_feasible_core_directions <-
    sum(final_target_feasible)
  model$build_params$association_flux_tolerance <- threshold
  model$build_params$target_flux_cell_type <- as.character(cell_type)
  model
}

.rc_corda_apply_target_epsilon <- function(
    model, epsilon, strict, cell_type) {
  .rc_corda2_apply_target_flux(
    model = model,
    flux_threshold = epsilon,
    strict = strict,
    cell_type = cell_type
  )
}
