#' Build medium-specific structural models and run directional LP scoring
#'
#' With `model_mode = "meta_module_gem"`, this stage is the only place where
#' FASTCORE is applied. For each medium scenario it constructs one union GEM
#' from the merged biological meta-module catalogue plus global FASTCORE support,
#' then reuses that exact model for every condition and metacell. Only union-GEM
#' construction accepts `model_params$completion_time_limit`; scoring LPs have
#' no time-limit parameter.
#'
#' @export
rc_regcompass_step_layer2 <- function(
    layer1, meta_modules, gem, medium_scenarios, outdir,
    model_mode = c("meta_module_gem", "full_gem"),
    layer2_args = list(), parallel = TRUE, BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  monitor <- .rc_step_monitor_start("layer2", outdir, progress)
  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)
  model_mode <- match.arg(model_mode)
  .rc_require_stage_class(
    meta_modules, "regcompass_meta_module_step", "meta_modules",
    "rc_regcompass_step_meta_modules"
  )
  if (!is.list(layer2_args)) {
    stop("`layer2_args` must be a list.", call. = FALSE)
  }
  allowed <- c(
    "model_params", "omega", "target_direction", "solver", "flux_threshold"
  )
  unknown <- setdiff(names(layer2_args), allowed)
  if (length(unknown)) {
    stop(
      "Unsupported `layer2_args`: ", paste(unknown, collapse = ", "),
      ". Scoring `time_limit`, sample aggregation, and sample columns have ",
      "been removed. Use only `layer2_args$model_params$completion_time_limit` ",
      "to limit global FASTCORE union-GEM construction.",
      call. = FALSE
    )
  }
  params <- meta_modules$workflow_params
  .rc_require_stage_gem(meta_modules, gem, "meta_modules")
  .rc_validate_layer1_stage(
    layer1, workflow_params = params, gem = gem, argument = "layer1"
  )
  medium_scenarios <- .rc_validate_shared_medium(medium_scenarios)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  layer2_args$model_params <- layer2_args$model_params %||% list()
  if (!is.list(layer2_args$model_params)) {
    stop("`layer2_args$model_params` must be a list.", call. = FALSE)
  }
  allowed_model_params <- c(
    "completion_time_limit", "fastcore_epsilon",
    "max_support_reactions", "strict"
  )
  unknown_model_params <- setdiff(
    names(layer2_args$model_params), allowed_model_params
  )
  if (length(unknown_model_params)) {
    stop(
      "Unsupported `layer2_args$model_params`: ",
      paste(unknown_model_params, collapse = ", "),
      ". Only union-GEM construction controls are accepted; ",
      "`time_limit` is not a scoring or model parameter.",
      call. = FALSE
    )
  }
  layer2_args$model_params$cache_dir <- file.path(
    outdir, "model_cache", model_mode
  )
  reserved <- intersect(names(layer2_args), c(
    "layer1", "gem", "mode", "unit", "reaction_membership",
    "core_reactions", "target_reactions", "medium_scenarios",
    "sample_col", "condition_col", "celltype_col", "BPPARAM",
    "parallel", "penalty_weights"
  ))
  if (length(reserved)) {
    stop(
      "`layer2_args` cannot override workflow fields: ",
      paste(reserved, collapse = ", "),
      call. = FALSE
    )
  }
  solver <- match.arg(
    as.character(layer2_args$solver %||% "highs"),
    c("highs", "gurobi", "glpk")
  )
  .rc_require_lp_solver(solver)
  catalogue <- meta_modules$merged_modules
  if (!is.list(catalogue) ||
      !is.data.frame(catalogue$merged_core_reactions) ||
      !is.data.frame(catalogue$merged_reaction_membership)) {
    stop(
      "The merged biological meta-module catalogue is incomplete.",
      call. = FALSE
    )
  }
  targets <- unique(as.character(
    catalogue$merged_core_reactions$reaction_id
  ))
  missing_expression <- setdiff(
    targets, rownames(layer1$reaction_expression)
  )
  if (length(missing_expression)) {
    stop(
      "Merged core reactions are absent from Layer 1 expression: ",
      paste(utils::head(missing_expression, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  defaults <- list(
    layer1 = layer1,
    gem = gem,
    target_reactions = targets,
    medium_scenarios = medium_scenarios,
    mode = model_mode,
    reaction_membership = if (identical(model_mode, "meta_module_gem")) {
      catalogue$merged_reaction_membership
    } else {
      NULL
    },
    core_reactions = if (identical(model_mode, "meta_module_gem")) {
      catalogue$merged_core_reactions
    } else {
      NULL
    },
    unit = "metacell",
    condition_col = params$condition_col,
    celltype_col = params$celltype_col,
    parallel = parallel,
    BPPARAM = BPPARAM
  )
  defaults[names(layer2_args)] <- NULL
  answer <- withCallingHandlers(
    do.call(rc_run_microcompass, c(defaults, layer2_args)),
    warning = function(w) {
      if (grepl(
        "Metacell-level scores are descriptive pseudo-observations",
        conditionMessage(w), fixed = TRUE
      )) invokeRestart("muffleWarning")
    }
  )
  answer$workflow_params <- params
  answer$gem_fingerprint <- .rc_stage_gem_fingerprint(gem)
  answer$source_core_reactions <- catalogue$merged_core_reactions
  answer$source_merged_reaction_membership <-
    catalogue$merged_reaction_membership
  answer$union_gem_policy <- if (identical(model_mode, "meta_module_gem")) {
    "one medium-specific union GEM; single global FASTCORE completion"
  } else {
    "shared full GEM; no union-GEM reconstruction"
  }
  class(answer) <- c("regcompass_layer2_step", "list")
  .rc_validate_layer2_stage(
    answer,
    layer1 = layer1,
    workflow_params = params,
    gem = gem,
    argument = "layer2"
  )
  answer <- .rc_step_monitor_finish(answer, monitor)
  rc_export_microcompass(answer, outdir)
  saveRDS(answer, file.path(outdir, "step_layer2.rds"))
  answer
}
