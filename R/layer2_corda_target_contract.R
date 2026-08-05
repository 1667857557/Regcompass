# Apply the Python CORDA2 target flux to closure and downstream scoring.

.rc_complete_celltype_medium_corda_gem_output_base <-
  .rc_complete_celltype_medium_corda_gem

.rc_corda2_apply_target_flux <- function(
    model, target_flux, strict, cell_type) {
  diagnostics <- as.data.frame(model$closure_diagnostics)
  required <- c(
    "reaction_id", "target_direction", "feasible", "vmax",
    "final_feasible", "final_vmax"
  )
  if (!all(required %in% colnames(diagnostics))) {
    stop("CORDA2 closure diagnostics cannot apply target flux.",
         call. = FALSE)
  }
  target_flux <- as.numeric(target_flux)
  if (length(target_flux) != 1L || !is.finite(target_flux) ||
      target_flux <= 0) {
    stop("CORDA2 target flux must be one positive finite number.",
         call. = FALSE)
  }
  parent_target_feasible <- diagnostics$feasible %in% TRUE &
    is.finite(diagnostics$vmax) & diagnostics$vmax >= target_flux
  final_target_feasible <- diagnostics$final_feasible %in% TRUE &
    is.finite(diagnostics$final_vmax) & diagnostics$final_vmax >= target_flux
  failed <- parent_target_feasible & !final_target_feasible
  diagnostics$parent_corda2_target_feasible <- parent_target_feasible
  diagnostics$final_corda2_target_feasible <- final_target_feasible
  diagnostics$corda2_target_flux <- target_flux
  diagnostics$completion_status <- ifelse(
    !parent_target_feasible,
    "parent_below_corda2_target_flux",
    ifelse(
      final_target_feasible,
      "corda2_retained_at_target_flux",
      "corda2_unresolved_at_target_flux"
    )
  )
  if (isTRUE(strict) && any(failed)) {
    bad <- paste(
      paste(
        diagnostics$reaction_id[failed],
        diagnostics$target_direction[failed],
        sep = ":"
      ),
      collapse = ", "
    )
    stop(
      "CORDA2 failed to retain parent-feasible core directions at target ",
      "flux in `", cell_type, "`: ", bad,
      call. = FALSE
    )
  }
  targets <- diagnostics[
    final_target_feasible,
    c("reaction_id", "target_direction"),
    drop = FALSE
  ]
  rownames(targets) <- NULL
  model$closure_diagnostics <- diagnostics
  model$target_directions <- targets
  model$target_status <- if (any(failed)) {
    "structurally_infeasible_at_corda2_target_flux"
  } else if (!nrow(targets)) {
    "parent_below_corda2_target_flux"
  } else {
    "ok"
  }
  model$build_params$corda2_target_flux <- target_flux
  model$build_params$n_parent_corda2_target_feasible_core_directions <-
    sum(parent_target_feasible)
  model$build_params$n_final_corda2_target_feasible_core_directions <-
    sum(final_target_feasible)
  model$build_params$association_flux_tolerance <-
    model$build_params$corda_options$flux_tolerance
  # Legacy audit aliases are retained for development objects created before
  # the public option was renamed to `corda2`.
  model$build_params$corda_target_epsilon <- target_flux
  model$build_params$n_parent_corda_epsilon_feasible_core_directions <-
    sum(parent_target_feasible)
  model$build_params$n_final_corda_epsilon_feasible_core_directions <-
    sum(final_target_feasible)
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
  target_flux <-
    args$corda_options$target_flux %||%
    args$corda_options$epsilon
  model <- .rc_corda2_apply_target_flux(
    model,
    target_flux = target_flux,
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
