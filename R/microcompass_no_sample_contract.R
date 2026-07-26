# Canonical RegCompass scoring is metacell resolved. Retain the audited LP
# engine while removing the historical sample-column interface and placeholder
# sample diagnostics from the active runtime contract.

.rc_run_microcompass_sample_legacy_core <- rc_run_microcompass

rc_run_microcompass <- function(
    layer1, gem, target_reactions = NULL,
    medium_table = NULL, medium_scenarios = NULL,
    mode = c("full_gem", "meta_module_gem"),
    reaction_membership = NULL, core_reactions = NULL,
    unit = "metacell", condition_col = "condition",
    celltype_col = "cell_type", model_params = list(),
    omega = 0.95,
    target_direction = c("both", "forward", "reverse"),
    parallel = TRUE,
    solver = c("highs", "gurobi", "glpk"),
    flux_threshold = 1e-8,
    BPPARAM = NULL) {
  if (!identical(as.character(unit), "metacell")) {
    stop(
      paste(
        "RegCompass scoring supports only `unit = \"metacell\"`.",
        "The historical sample-by-cell-type aggregation mode has been removed."
      ),
      call. = FALSE
    )
  }
  answer <- .rc_run_microcompass_sample_legacy_core(
    layer1 = layer1,
    gem = gem,
    target_reactions = target_reactions,
    medium_table = medium_table,
    medium_scenarios = medium_scenarios,
    mode = mode,
    reaction_membership = reaction_membership,
    core_reactions = core_reactions,
    unit = "metacell",
    condition_col = condition_col,
    sample_col = NULL,
    celltype_col = celltype_col,
    model_params = model_params,
    omega = omega,
    target_direction = target_direction,
    parallel = parallel,
    solver = solver,
    flux_threshold = flux_threshold,
    BPPARAM = BPPARAM
  )
  if (is.data.frame(answer$lp_diagnostics) &&
      "sample_id" %in% colnames(answer$lp_diagnostics)) {
    answer$lp_diagnostics$sample_id <- NULL
  }
  answer$params$unit <- "metacell"
  answer$params$aggregation <- "none"
  answer$params$sample_column <- NULL
  answer
}
