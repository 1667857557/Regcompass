# Exact CORDA2 Layer 2 public dispatch and completion metadata.

rc_regcompass_step_layer2 <- function(
    layer1, meta_modules, gem, medium_scenarios, outdir,
    model_mode = c("meta_module_gem", "full_gem"),
    layer2_args = list(), parallel = TRUE, BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  model_mode <- match.arg(model_mode)
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
    "model_completion",
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
  on.exit({
    rm(list = ls(.rc_layer2_completion_context, all.names = TRUE),
       envir = .rc_layer2_completion_context)
    list2env(previous, envir = .rc_layer2_completion_context)
  }, add = TRUE)
  answer <- .rc_regcompass_step_layer2_completion_base(
    layer1 = layer1,
    meta_modules = meta_modules,
    gem = gem,
    medium_scenarios = medium_scenarios,
    outdir = outdir,
    model_mode = model_mode,
    layer2_args = layer2_args,
    parallel = parallel,
    BPPARAM = BPPARAM,
    progress = progress
  )
  answer$params$model_completion <- if (isTRUE(is_corda2)) {
    "corda2"
  } else {
    "fastcore"
  }
  if (isTRUE(is_corda2)) {
    solver <- as.character(layer2_args$solver %||% "highs")
    feasibility_tolerance <-
      .rc_corda2_solver_feasibility_tolerance(solver)
    answer$completion_contract <- list(
      model_completion = "corda2",
      default_unchanged = FALSE,
      algorithm = corda_options$algorithm,
      source_fidelity = "exact_for_met_prod_NULL",
      python_reference = list(
        repository = "resendislab/corda",
        commit = corda_options$python_reference_commit,
        class = "CORDA2"
      ),
      confidence_levels = c(
        absent = -1L, unknown = 0L, low = 1L,
        medium = 2L, high = 3L
      ),
      confidence_mapping = paste(
        "high=merged core; medium=non-core module; low=optional",
        "high-evidence outside-module; absent=low evidence; unknown=remaining"
      ),
      constructor_controls = list(
        n = corda_options$redundancies,
        penalty_factor = corda_options$penalty_factor,
        support = corda_options$support
      ),
      fixed_source_constants = list(
        UPPER = 1e6,
        CI = 1.01,
        tflux = 1
      ),
      solver_configuration = list(
        solver = solver,
        feasibility_tolerance = feasibility_tolerance
      ),
      stage_update_policy = "python_serial_mutation_order",
      target_parallelism = FALSE,
      native_solver_acceleration = paste(
        "one persistent HiGHS C++ solver per CORDA2 instance with complete",
        "objective/bound updates and simplex-basis reuse"
      ),
      intentional_corrections = character(),
      options = corda_options
    )
    answer$params$structural_completion <- "corda2"
    answer$params$structural_completion_algorithm <- corda_options$algorithm
    answer$params$corda2_reference_commit <-
      corda_options$python_reference_commit
    answer$params$corda2_redundancies <- corda_options$redundancies
    answer$params$corda2_support <- corda_options$support
    answer$params$corda2_penalty_factor <- corda_options$penalty_factor
    answer$params$corda2_cost_increase <- 1.01
    answer$params$corda2_target_flux <- 1
    answer$params$corda2_feasibility_tolerance <- feasibility_tolerance
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
  } else {
    answer$completion_contract <- list(
      model_completion = "fastcore",
      default_unchanged = TRUE,
      algorithm = "add_only_compact_FASTCORE"
    )
  }
  rc_export_microcompass(answer, outdir)
  model_dir <- file.path(outdir, "03_models")
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    answer$completion_contract,
    file.path(model_dir, "model_completion_contract.rds")
  )
  saveRDS(answer, file.path(outdir, "step_layer2.rds"))
  answer
}
