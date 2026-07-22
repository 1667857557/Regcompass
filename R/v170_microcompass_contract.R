# Canonical microCOMPASS entry: the explicit `unit` argument is authoritative.
.rc_run_microcompass_with_explicit_unit <- rc_run_microcompass
.rc_run_microcompass_v170 <- function(
    layer1, gem, target_reactions = NULL,
    medium_table = NULL, medium_scenarios = NULL,
    mode = c("full_gem", "meta_module_gem"),
    reaction_membership = NULL, core_reactions = NULL,
    unit = c("sample_celltype", "metacell"),
    condition_col = "condition", sample_col = "sample_id",
    celltype_col = "cell_type", model_params = list(),
    omega = 0.95,
    target_direction = c("both", "forward", "reverse"),
    parallel = TRUE,
    solver = c("highs", "gurobi", "glpk"),
    time_limit = 60, flux_threshold = 1e-8,
    BPPARAM = NULL) {
  solver <- match.arg(solver)
  .rc_require_lp_solver(solver)
  hidden_unit <- getOption("RegCompassR.inference_unit", NULL)
  if (!is.null(hidden_unit)) {
    warning(
      "The retired `RegCompassR.inference_unit` option is ignored; use the explicit `unit` argument.",
      call. = FALSE
    )
  }
  old_option <- options("RegCompassR.inference_unit")
  options(RegCompassR.inference_unit = NULL)
  on.exit(options(old_option), add = TRUE)
  .rc_run_microcompass_with_explicit_unit(
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
    omega = omega,
    target_direction = target_direction,
    parallel = parallel,
    solver = solver,
    time_limit = time_limit,
    flux_threshold = flux_threshold,
    BPPARAM = BPPARAM
  )
}
rc_run_microcompass <- .rc_run_microcompass_v170
