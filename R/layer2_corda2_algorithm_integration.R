# Attach exact Python CORDA2 reconstruction metadata to Layer 2 models.

.rc_complete_celltype_medium_corda_gem_before_corda2 <-
  .rc_complete_celltype_medium_corda_gem

.rc_complete_celltype_medium_corda_gem <- function(...) {
  model <- .rc_complete_celltype_medium_corda_gem_before_corda2(...)
  reconstruction <- model$corda_reconstruction
  model$build_params$strategy <- "celltype_medium_python_corda2_exact"
  model$build_params$algorithm <- reconstruction$algorithm
  model$build_params$completion_stage <-
    "python_CORDA2_exact_after_confidence_mapping"
  model$build_params$stage_update_policy <- reconstruction$stage_update_policy
  model$build_params$python_reference_commit <-
    reconstruction$python_reference_commit
  model$build_params$python_source_semantics <- reconstruction$source_semantics
  model$build_params$corda2_redundancies <-
    model$build_params$corda_options$redundancies
  model$build_params$corda2_support <-
    model$build_params$corda_options$support
  model$build_params$corda2_penalty_factor <-
    model$build_params$corda_options$penalty_factor
  model$build_params$corda2_cost_increase <- 1.01
  model$build_params$corda2_target_flux <- 1
  model$build_params$corda2_upper_bound <- 1e6
  model$build_params$corda2_feasibility_tolerance <-
    model$build_params$corda_options$feasibility_tolerance
  model$corda2_contract <- list(
    implementation = "exact resendislab/corda Python CORDA2 semantics",
    supported_scope = "met_prod = NULL",
    reference_repository = "resendislab/corda",
    reference_commit = reconstruction$python_reference_commit,
    confidence_levels = c(absent = -1L, unknown = 0L,
                          low = 1L, medium = 2L, high = 3L),
    redundant_path_cost_increase = 1.01,
    absent_penalty_factor = model$build_params$corda2_penalty_factor,
    support_threshold = model$build_params$corda2_support,
    maximum_redundant_paths = model$build_params$corda2_redundancies,
    target_flux = 1,
    upper_bound = 1e6,
    feasibility_tolerance =
      model$build_params$corda2_feasibility_tolerance,
    source_semantics = reconstruction$source_semantics
  )
  model
}

.rc_complete_celltype_medium_corda_like_gem <-
  .rc_complete_celltype_medium_corda_gem
