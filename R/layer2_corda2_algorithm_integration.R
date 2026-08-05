# Attach exact Python CORDA2 reconstruction metadata to Layer 2 models.

.rc_complete_celltype_medium_corda_gem_before_corda2 <-
  .rc_complete_celltype_medium_corda_gem

.rc_complete_celltype_medium_corda_gem <- function(...) {
  model <- .rc_complete_celltype_medium_corda_gem_before_corda2(...)
  reconstruction <- model$corda_reconstruction
  options <- model$build_params$corda_options
  model$build_params$strategy <- "celltype_medium_python_corda2_exact"
  model$build_params$algorithm <- reconstruction$algorithm
  model$build_params$completion_stage <-
    "python_CORDA2_exact_after_confidence_mapping"
  model$build_params$stage_update_policy <- reconstruction$stage_update_policy
  model$build_params$python_reference_commit <-
    reconstruction$python_reference_commit
  model$build_params$python_source_semantics <- reconstruction$source_semantics

  # Canonical names reproduce CORDA.__init__(model, confidence, met_prod=None,
  # n=3, penalty_factor=100, support=5). Legacy fields remain read-only output
  # aliases for existing RegCompass result consumers.
  model$build_params$corda2_args <- list(
    met_prod = options$met_prod,
    n = options$n,
    penalty_factor = options$penalty_factor,
    support = options$support
  )
  model$build_params$corda2_met_prod <- options$met_prod
  model$build_params$corda2_n <- options$n
  model$build_params$corda2_penalty_factor <- options$penalty_factor
  model$build_params$corda2_support <- options$support
  model$build_params$corda2_redundancies <- options$n
  model$build_params$corda2_cost_increase <- 1.01
  model$build_params$corda2_target_flux <- 1
  model$build_params$corda2_upper_bound <- 1e6
  model$build_params$corda2_feasibility_tolerance <-
    options$feasibility_tolerance
  model$build_params$corda2_solver_time_limit <-
    reconstruction$solver_time_limit %||% Inf
  model$build_params$regcompass_time_limit_ignored_for_corda2 <-
    reconstruction$requested_regcompass_time_limit_ignored %||% NA_real_

  model$corda2_contract <- list(
    implementation = "exact resendislab/corda Python CORDA2 semantics",
    supported_scope = "met_prod = NULL",
    reference_repository = "resendislab/corda",
    reference_commit = reconstruction$python_reference_commit,
    constructor_signature = c(
      "model", "confidence", "met_prod", "n", "penalty_factor", "support"
    ),
    constructor_args = model$build_params$corda2_args,
    confidence_levels = c(absent = -1L, unknown = 0L,
                          low = 1L, medium = 2L, high = 3L),
    fixed_constants = c(CI = 1.01, tflux = 1, UPPER = 1e6),
    solver_time_limit = reconstruction$solver_time_limit %||% Inf,
    feasibility_tolerance = options$feasibility_tolerance,
    source_semantics = reconstruction$source_semantics
  )
  model
}

.rc_complete_celltype_medium_corda_like_gem <-
  .rc_complete_celltype_medium_corda_gem
