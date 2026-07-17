# Preserve public validation order and named diagnostics after the required
# result-level corrections.

.rc_required_result_run_microcompass <- rc_run_microcompass

rc_run_microcompass <- function(
    layer1, gem, target_reactions = NULL,
    medium_table = NULL, medium_scenarios = NULL,
    mode = c("full_gem", "meta_module_gem"),
    reaction_membership = NULL, core_reactions = NULL,
    unit = c("sample_celltype", "metacell"),
    condition_col = "condition", sample_col = "sample_id",
    celltype_col = "cell_type", model_params = list(),
    penalty_weights = c(expr = 1.0, confidence = 0.5, missing = 1.0),
    omega = 0.95,
    target_direction = c("both", "forward", "reverse"),
    parallel = TRUE,
    solver = c("highs", "gurobi", "glpk"),
    time_limit = 60, flux_threshold = 1e-8,
    BPPARAM = NULL) {
  # Invalid public arguments must fail before attempting to inspect the GEM.
  mode <- match.arg(mode)
  unit <- match.arg(unit)
  target_direction <- match.arg(target_direction)
  solver <- match.arg(solver)

  .rc_required_result_run_microcompass(
    layer1 = layer1,
    gem = gem,
    target_reactions = target_reactions,
    medium_table = medium_table,
    medium_scenarios = medium_scenarios,
    mode = mode,
    reaction_membership = reaction_membership,
    core_reactions = core_reactions,
    unit = unit,
    condition_col = condition_col,
    sample_col = sample_col,
    celltype_col = celltype_col,
    model_params = model_params,
    penalty_weights = penalty_weights,
    omega = omega,
    target_direction = target_direction,
    parallel = parallel,
    solver = solver,
    time_limit = time_limit,
    flux_threshold = flux_threshold,
    BPPARAM = BPPARAM
  )
}

rc_compute_multiome_penalty <- function(...) {
  answer <- .rc_compute_multiome_penalty_core(...)
  reaction_ids <- rownames(answer$penalty)
  named_components <- c(
    "gpr_missing_fraction", "role", "role_source", "role_confidence",
    "role_override_flag", "support_penalty_used"
  )
  for (name in named_components) {
    value <- answer$components[[name]]
    if (!is.null(value) && is.null(dim(value)) &&
        length(value) == length(reaction_ids)) {
      names(value) <- reaction_ids
      answer$components[[name]] <- value
    }
  }
  answer$evidence_policy <- "penalty_only"
  answer$evidence_description <- paste(
    "This is not the original COMPASS expression-neighbourhood penalty.",
    "RegCompass uses a COMPASS-like inverse-support expression term;",
    "multiome evidence modifies the LP objective penalty only and does not",
    "directly change stoichiometry or internal reaction bounds."
  )
  answer
}
