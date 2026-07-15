.rc_build_meta_module_gem_before_status_fix <- rc_build_meta_module_gem

rc_build_meta_module_gem <- function(...) {
  model <- .rc_build_meta_module_gem_before_status_fix(...)
  diagnostics <- model$closure_diagnostics
  if (!is.data.frame(diagnostics) || !nrow(diagnostics) ||
      !all(c("target_direction", "completion_status") %in% colnames(diagnostics))) {
    return(model)
  }

  parent_blocked <- diagnostics$target_direction %in% c("forward", "reverse") &
    diagnostics$completion_status == "no_allowed_direction"
  diagnostics$completion_status[parent_blocked] <- "parent_blocked"
  model$closure_diagnostics <- diagnostics

  failed <- diagnostics$feasible %in% TRUE &
    !(diagnostics$final_feasible %in% TRUE)
  n_no_direction <- sum(diagnostics$completion_status == "no_allowed_direction")
  n_parent_blocked <- sum(diagnostics$completion_status == "parent_blocked")
  model$target_status <- if (any(failed)) {
    "structurally_infeasible"
  } else if (nrow(diagnostics) > 0L && n_no_direction == nrow(diagnostics)) {
    "no_allowed_direction"
  } else if (n_no_direction > 0L) {
    "partial_no_allowed_direction"
  } else if (nrow(diagnostics) > 0L && n_parent_blocked == nrow(diagnostics)) {
    "parent_blocked"
  } else if (n_parent_blocked > 0L) {
    "partial_parent_blocked"
  } else {
    "ok"
  }
  model
}
