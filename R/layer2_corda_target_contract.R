# Prepare post-CORDA2 core targets without any post-build LP solve.

.rc_corda2_prepare_scoring_targets <- function(
    model, core_reactions,
    target_direction = c("both", "forward", "reverse"),
    strict = TRUE, cell_type = NA_character_) {
  target_direction <- match.arg(target_direction)
  validated <- rc_validate_gem(model)
  core <- unique(trimws(as.character(core_reactions)))
  core <- core[!is.na(core) & nzchar(core)]
  if (!length(core)) {
    stop("CORDA2 scoring requires at least one core reaction.", call. = FALSE)
  }

  missing_core <- setdiff(core, validated$reactions)
  if (length(missing_core)) {
    stop(
      "CORDA2 final GEM is missing required core reactions: ",
      paste(utils::head(missing_core, 10L), collapse = ", "),
      ". Core reactions are an immutable structural backbone.",
      call. = FALSE
    )
  }

  diagnostics <- rc_prepare_directional_targets(
    model,
    core,
    target_direction = target_direction
  )
  targets <- diagnostics[
    diagnostics$target_direction %in% c("forward", "reverse"),
    c("reaction_id", "target_direction"),
    drop = FALSE
  ]
  rownames(targets) <- NULL

  model$required_core_reactions <- core
  model$core_direction_diagnostics <- diagnostics
  model$target_directions <- targets
  model$target_status <- if (!nrow(targets)) {
    "no_bound_allowed_core_direction"
  } else {
    "ready_for_compass_vmax"
  }

  if (!is.list(model$build_params)) model$build_params <- list()
  model$build_params$strict_requested <- strict
  model$build_params$strict_used_for_reconstruction <- FALSE
  model$build_params$n_required_core_reactions <- length(core)
  model$build_params$n_retained_core_reactions <- sum(
    core %in% validated$reactions
  )
  model$build_params$n_core_scoring_directions <- nrow(targets)
  model$build_params$post_reconstruction_closure_lp <- FALSE
  model$build_params$target_feasibility <- paste(
    "no post-CORDA2 closure LP; microCOMPASS computes directional vmax once",
    "on the reconstructed GEM and skips the penalty LP when infeasible"
  )
  model$build_params$target_flux_cell_type <- as.character(cell_type)
  model
}
