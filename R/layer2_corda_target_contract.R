# Apply fixed Python CORDA2 tflux=1 to RegCompass post-build scoring targets.

.rc_complete_celltype_medium_corda_gem_output_base <-
  .rc_complete_celltype_medium_corda_gem

.rc_corda2_apply_target_flux <- function(
    model, target_flux = 1, strict, cell_type) {
  diagnostics <- as.data.frame(model$closure_diagnostics)
  required <- c(
    "reaction_id", "target_direction", "feasible", "vmax",
    "final_feasible", "final_vmax"
  )
  if (!all(required %in% colnames(diagnostics))) {
    stop("CORDA2 closure diagnostics cannot apply tflux.", call. = FALSE)
  }
  if (!identical(as.numeric(target_flux), 1)) {
    stop("Exact Python CORDA2 fixes `tflux` at 1.", call. = FALSE)
  }
  parent_target_feasible <- diagnostics$feasible %in% TRUE &
    is.finite(diagnostics$vmax) & diagnostics$vmax >= 1
  final_target_feasible <- diagnostics$final_feasible %in% TRUE &
    is.finite(diagnostics$final_vmax) & diagnostics$final_vmax >= 1
  failed <- parent_target_feasible & !final_target_feasible
  diagnostics$parent_corda2_tflux_feasible <- parent_target_feasible
  diagnostics$final_corda2_tflux_feasible <- final_target_feasible
  diagnostics$corda2_tflux <- 1
  diagnostics$completion_status <- ifelse(
    !parent_target_feasible,
    "parent_below_corda2_tflux",
    ifelse(
      final_target_feasible,
      "corda2_retained_at_tflux",
      "corda2_unresolved_at_tflux"
    )
  )
  # Exact Python CORDA2 does not stop when a target is impossible; it records
  # the target in `impossible` and completes the build. Do not elevate this
  # RegCompass post-build closure diagnostic into an algorithmic error.
  targets <- diagnostics[
    final_target_feasible,
    c("reaction_id", "target_direction"),
    drop = FALSE
  ]
  rownames(targets) <- NULL
  model$closure_diagnostics <- diagnostics
  model$target_directions <- targets
  model$target_status <- if (any(failed)) {
    "structurally_infeasible_at_corda2_tflux"
  } else if (!nrow(targets)) {
    "parent_below_corda2_tflux"
  } else {
    "ok"
  }
  model$build_params$corda2_target_flux <- 1
  model$build_params$strict_requested <- strict
  model$build_params$strict_used_for_reconstruction <- FALSE
  model$build_params$n_parent_corda2_tflux_feasible_core_directions <-
    sum(parent_target_feasible)
  model$build_params$n_final_corda2_tflux_feasible_core_directions <-
    sum(final_target_feasible)
  model$build_params$association_flux_tolerance <-
    model$build_params$corda_options$feasibility_tolerance
  model
}

.rc_corda_apply_target_epsilon <- function(
    model, epsilon, strict, cell_type) {
  .rc_corda2_apply_target_flux(
    model = model,
    target_flux = epsilon,
    strict = strict,
    cell_type = cell_type
  )
}

.rc_complete_celltype_medium_corda_gem <- function(...) {
  args <- list(...)
  model <- do.call(.rc_complete_celltype_medium_corda_gem_output_base, args)
  model <- .rc_corda2_apply_target_flux(
    model,
    target_flux = args$corda_options$target_flux,
    strict = args$strict,
    cell_type = as.character(args$cell_type)
  )
  .rc_validate_corda_union_model(
    model,
    cell_type = as.character(args$cell_type)
  )
  model
}

.rc_complete_celltype_medium_corda_like_gem <-
  .rc_complete_celltype_medium_corda_gem
