# Apply the fixed original CORDA2 flux threshold to post-build scoring targets.

.rc_corda2_apply_target_flux <- function(
    model, flux_threshold = 1e-7, strict, cell_type) {
  diagnostics <- as.data.frame(model$closure_diagnostics)
  required <- c("reaction_id", "target_direction", "feasible", "vmax")
  if (!all(required %in% colnames(diagnostics))) {
    stop("CORDA2 closure diagnostics are incomplete.", call. = FALSE)
  }
  threshold <- as.numeric(flux_threshold)
  if (length(threshold) != 1L || !is.finite(threshold) || threshold <= 0) {
    stop("CORDA2 flux threshold must be positive.", call. = FALSE)
  }
  target_feasible <- diagnostics$feasible %in% TRUE &
    is.finite(diagnostics$vmax) & diagnostics$vmax >= threshold
  diagnostics$corda2_feasible <- target_feasible
  diagnostics$corda2_flux_threshold <- threshold
  diagnostics$completion_status <- ifelse(
    target_feasible, "corda2_retained", "corda2_blocked"
  )
  targets <- diagnostics[
    target_feasible,
    c("reaction_id", "target_direction"),
    drop = FALSE
  ]
  rownames(targets) <- NULL
  model$closure_diagnostics <- diagnostics
  model$target_directions <- targets
  model$target_status <- if (!nrow(targets)) {
    "no_feasible_core_direction"
  } else {
    "ok"
  }
  model$build_params$strict_requested <- strict
  model$build_params$strict_used_for_reconstruction <- FALSE
  model$build_params$n_corda2_tested_core_directions <- nrow(diagnostics)
  model$build_params$n_corda2_feasible_core_directions <- sum(target_feasible)
  model$build_params$association_flux_tolerance <- threshold
  model$build_params$target_flux_cell_type <- as.character(cell_type)
  model
}
