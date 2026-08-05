# Final public/runtime metadata contract for exact Python CORDA2 semantics.

.rc_is_corda2_options <- function(options) {
  is.list(options) && identical(
    as.character(options$algorithm %||% ""),
    "resendislab_python_CORDA2_c02e06d_exact_semantics"
  )
}

.rc_build_celltype_medium_union_gem_cache_exact_base <-
  .rc_build_celltype_medium_union_gem_cache

.rc_build_celltype_medium_union_gem_cache <- function(...) {
  cache <- .rc_build_celltype_medium_union_gem_cache_exact_base(...)
  if (!identical(attr(cache, "completion_method"), "corda2")) return(cache)
  for (name in names(cache)) {
    cache[[name]]$build_strategy <- "celltype_medium_python_corda2_exact"
  }
  summary <- attr(cache, "summary")
  if (is.data.frame(summary) && nrow(summary)) {
    summary$build_strategy <- "celltype_medium_python_corda2_exact"
    summary$completion_stage <- "python_CORDA2_exact_after_confidence_mapping"
    attr(cache, "summary") <- summary
  }
  attr(cache, "structural_parallel_task") <- if (
    identical(attr(cache, "structural_parallel_task"),
              "cell_type_x_medium_outer_parallel")
  ) {
    "cell_type_x_medium_outer_parallel_python_instances"
  } else {
    "serial_within_each_python_corda2_instance"
  }
  cache
}

.rc_regcompass_step_layer2_exact_public_base <- rc_regcompass_step_layer2

rc_regcompass_step_layer2 <- function(
    layer1, meta_modules, gem, medium_scenarios, outdir,
    model_mode = c("meta_module_gem", "full_gem"),
    layer2_args = list(), parallel = TRUE, BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  answer <- .rc_regcompass_step_layer2_exact_public_base(
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
  model_params <- if (is.list(layer2_args)) {
    layer2_args$model_params %||% list()
  } else {
    list()
  }
  requested <- as.character(model_params$model_completion %||% "fastcore")
  if (!requested %in% c("corda2", "corda", "corda_like")) return(answer)
  options <- .rc_layer2_corda_options(model_params)
  constructor_args <- list(
    met_prod = options$met_prod,
    n = options$n,
    penalty_factor = options$penalty_factor,
    support = options$support
  )
  answer$params$model_completion <- "corda2"
  answer$params$structural_completion <- "corda2"
  answer$params$structural_completion_algorithm <- options$algorithm
  answer$params$corda2_reference_commit <- options$python_reference_commit
  answer$params$corda2_args <- constructor_args
  answer$params$corda2_met_prod <- options$met_prod
  answer$params$corda2_n <- options$n
  answer$params$corda2_penalty_factor <- options$penalty_factor
  answer$params$corda2_support <- options$support
  # Compatibility output alias; canonical Python name is `n`.
  answer$params$corda2_redundancies <- options$n
  answer$params$corda2_cost_increase <- 1.01
  answer$params$corda2_target_flux <- 1
  solver <- as.character(layer2_args$solver %||% "highs")
  feasibility_tolerance <-
    .rc_corda2_solver_feasibility_tolerance(solver)
  answer$params$corda2_feasibility_tolerance <- feasibility_tolerance
  answer$params$corda2_solver_time_limit <- Inf
  answer$completion_contract <- list(
    model_completion = "corda2",
    default_unchanged = FALSE,
    algorithm = options$algorithm,
    source_fidelity = "exact_for_met_prod_NULL",
    python_reference = list(
      repository = "resendislab/corda",
      commit = options$python_reference_commit,
      class = "CORDA"
    ),
    constructor_signature = c(
      "model", "confidence", "met_prod", "n", "penalty_factor", "support"
    ),
    constructor_args = constructor_args,
    fixed_source_constants = list(
      CI = 1.01,
      tflux = 1,
      UPPER = 1e6
    ),
    confidence_levels = c(
      absent = -1L, unknown = 0L, low = 1L,
      medium = 2L, high = 3L
    ),
    solver_configuration = list(
      solver = solver,
      threads = 1L,
      algorithm = if (identical(solver, "highs")) "simplex" else NA_character_,
      feasibility_tolerance = feasibility_tolerance,
      time_limit = Inf,
      option_round_trip_required = identical(solver, "highs")
    ),
    exact_source_behaviors = c(
      "ordered COBRA reaction-bound setters during construction",
      "both direction variables exist for every reaction",
      "target assessment does not block the opposite direction",
      "both direction costs use the forward variable confidence",
      "targets and stage-2 medium variables are processed serially",
      "stage-2 medium test minimizes a positive target coefficient",
      "no reconstruction time limit is added"
    ),
    intentional_corrections = character(),
    supported_scope = "met_prod = NULL"
  )
  answer$union_gem_policy <- paste(
    "one exact Python-CORDA2 instance per cell type and medium;",
    "parallelism is only across independent model instances, while each",
    "instance preserves Python target and variable mutation order"
  )
  answer$method <- paste(
    "microCOMPASS directional LP on cell-type-specific medium models",
    "reconstructed with exact resendislab Python CORDA2 semantics"
  )
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
