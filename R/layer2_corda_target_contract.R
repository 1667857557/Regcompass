# Apply the paper epsilon to CORDA target feasibility and downstream scoring.

.rc_complete_celltype_medium_corda_gem_output_base <-
  .rc_complete_celltype_medium_corda_gem

.rc_corda_apply_target_epsilon <- function(model, epsilon, strict, cell_type) {
  diagnostics <- as.data.frame(model$closure_diagnostics)
  required <- c(
    "reaction_id", "target_direction", "feasible", "vmax",
    "final_feasible", "final_vmax"
  )
  if (!all(required %in% colnames(diagnostics))) {
    stop("CORDA closure diagnostics cannot apply the paper epsilon.",
         call. = FALSE)
  }
  parent_epsilon_feasible <- diagnostics$feasible %in% TRUE &
    is.finite(diagnostics$vmax) & diagnostics$vmax >= epsilon
  final_epsilon_feasible <- diagnostics$final_feasible %in% TRUE &
    is.finite(diagnostics$final_vmax) & diagnostics$final_vmax >= epsilon
  failed <- parent_epsilon_feasible & !final_epsilon_feasible
  diagnostics$parent_corda_epsilon_feasible <- parent_epsilon_feasible
  diagnostics$final_corda_epsilon_feasible <- final_epsilon_feasible
  diagnostics$corda_epsilon <- epsilon
  diagnostics$completion_status <- ifelse(
    !parent_epsilon_feasible,
    "parent_below_corda_epsilon",
    ifelse(
      final_epsilon_feasible,
      "corda_retained_at_epsilon",
      "corda_unresolved_at_epsilon"
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
      "CORDA failed to retain parent-feasible HC directions at epsilon in `",
      cell_type, "`: ", bad,
      call. = FALSE
    )
  }
  targets <- diagnostics[
    final_epsilon_feasible,
    c("reaction_id", "target_direction"),
    drop = FALSE
  ]
  rownames(targets) <- NULL
  model$closure_diagnostics <- diagnostics
  model$target_directions <- targets
  model$target_status <- if (any(failed)) {
    "structurally_infeasible_at_corda_epsilon"
  } else if (!nrow(targets)) {
    "parent_below_corda_epsilon"
  } else {
    "ok"
  }
  model$build_params$corda_target_epsilon <- epsilon
  model$build_params$n_parent_corda_epsilon_feasible_core_directions <-
    sum(parent_epsilon_feasible)
  model$build_params$n_final_corda_epsilon_feasible_core_directions <-
    sum(final_epsilon_feasible)
  model$build_params$association_flux_tolerance <-
    model$build_params$corda_options$flux_tolerance
  model
}

.rc_complete_celltype_medium_corda_gem <- function(...) {
  args <- list(...)
  model <- do.call(.rc_complete_celltype_medium_corda_gem_output_base, args)
  model <- .rc_corda_apply_target_epsilon(
    model,
    epsilon = args$corda_options$epsilon,
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
