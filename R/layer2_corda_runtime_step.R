# CORDA2 preparation and completion helpers used directly by Layer 2.

.rc_layer2_prepare_completion <- function(
    layer1, meta_modules, model_mode, layer2_args) {
  if (!is.list(layer2_args)) {
    stop("`layer2_args` must be a list.", call. = FALSE)
  }
  model_params <- layer2_args$model_params %||% list()
  corda_options <- .rc_layer2_corda_options(model_params)
  is_corda2 <- .rc_is_corda2_options(corda_options)
  if (isTRUE(is_corda2) && !identical(model_mode, "meta_module_gem")) {
    stop(
      "`model_completion = \"corda2\"` is available only with ",
      "`model_mode = \"meta_module_gem\"`.",
      call. = FALSE
    )
  }

  extracted <- c(
    "model_completion", "corda2_args",
    "corda2_penalty_factor", "corda_penalty_factor", "corda_gamma",
    "corda2_cost_increase", "corda_cost_increase", "corda_kappa",
    "corda2_target_flux", "corda_tflux", "corda_epsilon",
    "corda2_redundancies", "corda_n",
    "corda2_support", "corda_support", "corda_p",
    "corda2_flux_tolerance", "corda_flux_tolerance", "corda_seed",
    "corda_medium_confidence_threshold",
    "corda_negative_confidence_threshold",
    "corda_regulatory_weight",
    "corda_include_evidence_outside_modules",
    "corda_max_medium_confidence_reactions",
    "corda_other_penalty", "corda_negative_penalty"
  )
  clean_params <- model_params[setdiff(names(model_params), extracted)]
  layer2_args$model_params <- clean_params

  previous <- as.list(.rc_layer2_completion_context)
  .rc_layer2_completion_context$active <- TRUE
  .rc_layer2_completion_context$model_completion <-
    if (isTRUE(is_corda2)) "corda2" else "fastcore"
  .rc_layer2_completion_context$corda_options <- corda_options
  .rc_layer2_completion_context$reaction_evidence <- if (isTRUE(is_corda2)) {
    .rc_layer2_corda_reaction_evidence(
      layer1,
      meta_modules,
      regulatory_weight = corda_options$regulatory_weight
    )
  } else {
    NULL
  }

  list(
    layer2_args = layer2_args,
    corda_options = corda_options,
    is_corda2 = is_corda2,
    previous_context = previous
  )
}

.rc_layer2_restore_completion <- function(previous_context) {
  rm(
    list = ls(.rc_layer2_completion_context, all.names = TRUE),
    envir = .rc_layer2_completion_context
  )
  list2env(previous_context, envir = .rc_layer2_completion_context)
  invisible(NULL)
}

.rc_layer2_finalize_completion <- function(
    answer, corda_options, is_corda2, solver) {
  answer$params$model_completion <- if (isTRUE(is_corda2)) {
    "corda2"
  } else {
    "fastcore"
  }
  if (!isTRUE(is_corda2)) {
    answer$completion_contract <- list(
      model_completion = "fastcore",
      default_unchanged = TRUE,
      algorithm = "add_only_compact_FASTCORE"
    )
    return(answer)
  }

  feasibility_tolerance <-
    .rc_corda2_solver_feasibility_tolerance(solver)
  constructor_args <- list(
    met_prod = corda_options$met_prod,
    n = corda_options$n,
    penalty_factor = corda_options$penalty_factor,
    support = corda_options$support
  )
  answer$completion_contract <- list(
    model_completion = "corda2",
    default_unchanged = FALSE,
    algorithm = corda_options$algorithm,
    source_fidelity = "exact_for_met_prod_NULL",
    python_reference = list(
      repository = "resendislab/corda",
      commit = corda_options$python_reference_commit,
      class = "CORDA"
    ),
    constructor_signature = c(
      "model", "confidence", "met_prod", "n", "penalty_factor", "support"
    ),
    constructor_args = constructor_args,
    confidence_levels = c(
      absent = -1L, unknown = 0L, low = 1L,
      medium = 2L, high = 3L
    ),
    fixed_source_constants = list(UPPER = 1e6, CI = 1.01, tflux = 1),
    solver_configuration = list(
      solver = solver,
      threads = 1L,
      algorithm = if (identical(solver, "highs")) "simplex" else NA_character_,
      feasibility_tolerance = feasibility_tolerance,
      time_limit = Inf,
      option_round_trip_required = identical(solver, "highs")
    ),
    stage_update_policy = "python_serial_mutation_order",
    target_parallelism = FALSE,
    supported_scope = "met_prod = NULL",
    intentional_corrections = character(),
    options = corda_options
  )
  answer$params$structural_completion <- "corda2"
  answer$params$structural_completion_algorithm <- corda_options$algorithm
  answer$params$corda2_reference_commit <-
    corda_options$python_reference_commit
  answer$params$corda2_args <- constructor_args
  answer$params$corda2_met_prod <- corda_options$met_prod
  answer$params$corda2_n <- corda_options$n
  answer$params$corda2_redundancies <- corda_options$n
  answer$params$corda2_support <- corda_options$support
  answer$params$corda2_penalty_factor <- corda_options$penalty_factor
  answer$params$corda2_cost_increase <- 1.01
  answer$params$corda2_target_flux <- 1
  answer$params$corda2_feasibility_tolerance <- feasibility_tolerance
  answer$params$corda2_solver_time_limit <- Inf
  answer$params$corda2_inner_target_parallelism <- FALSE
  answer$union_gem_policy <- paste(
    "one exact pinned Python-CORDA2 reconstruction per cell type and",
    "medium; each instance preserves serial target and variable mutation",
    "order; parallelism is limited to independent model instances"
  )
  answer$method <- paste(
    "microCOMPASS directional LP on cell-type-specific medium models",
    "reconstructed with exact pinned Python CORDA2 semantics"
  )
  answer
}
