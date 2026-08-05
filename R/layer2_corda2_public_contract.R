# Public Layer 2 contract for the corrected Python CORDA2 implementation.

.rc_regcompass_step_layer2_before_corda2_public <- rc_regcompass_step_layer2

rc_regcompass_step_layer2 <- function(
    layer1, meta_modules, gem, medium_scenarios, outdir,
    model_mode = c("meta_module_gem", "full_gem"),
    layer2_args = list(), parallel = TRUE, BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  model_mode <- match.arg(model_mode)
  if (!is.list(layer2_args)) {
    stop("`layer2_args` must be a list.", call. = FALSE)
  }
  original_params <- layer2_args$model_params %||% list()
  requested <- as.character(
    original_params$model_completion %||% "fastcore"
  )
  is_corda2 <- length(requested) == 1L && !is.na(requested) &&
    requested %in% c("corda2", "corda", "corda_like")
  translated <- original_params
  if (is_corda2) {
    translated$model_completion <- "corda"
    translated$corda_gamma <-
      original_params$corda2_penalty_factor %||%
      original_params$corda_penalty_factor %||%
      original_params$corda_gamma %||% 100
    translated$corda_kappa <-
      original_params$corda2_cost_increase %||%
      original_params$corda_cost_increase %||%
      original_params$corda_kappa %||% 1.01
    translated$corda_epsilon <-
      original_params$corda2_target_flux %||%
      original_params$corda_tflux %||%
      original_params$corda_epsilon %||% 1
    translated$corda_n <-
      original_params$corda2_redundancies %||%
      original_params$corda_n %||% 3L
    translated$corda_p <-
      original_params$corda2_support %||%
      original_params$corda_support %||%
      original_params$corda_p %||% 5L
    translated$corda_flux_tolerance <-
      original_params$corda2_flux_tolerance %||%
      original_params$corda_flux_tolerance %||% 1e-8
    translated[c(
      "corda2_penalty_factor", "corda_penalty_factor",
      "corda2_cost_increase", "corda_cost_increase",
      "corda2_target_flux", "corda_tflux",
      "corda2_redundancies", "corda2_support", "corda_support",
      "corda2_flux_tolerance"
    )] <- NULL
  }
  layer2_args$model_params <- translated
  answer <- .rc_regcompass_step_layer2_before_corda2_public(
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
  if (is_corda2) {
    options <- .rc_layer2_corda_options(original_params)
    answer$params$model_completion <- "corda2"
    answer$params$structural_completion <- "corda2"
    answer$params$structural_completion_algorithm <-
      "resendislab_python_CORDA2_corrected_redundant_path_assessment"
    answer$params$corda2_reference_commit <-
      options$python_reference_commit
    answer$params$corda2_redundancies <- options$redundancies
    answer$params$corda2_support <- options$support
    answer$params$corda2_penalty_factor <- options$penalty_factor
    answer$params$corda2_cost_increase <- options$cost_increase
    answer$params$corda2_target_flux <- options$target_flux
    answer$completion_contract$model_completion <- "corda2"
    answer$completion_contract$algorithm <-
      "resendislab_python_CORDA2_corrected_redundant_path_assessment"
    answer$completion_contract$python_reference <- list(
      repository = "resendislab/corda",
      commit = options$python_reference_commit,
      class = "CORDA2"
    )
    answer$completion_contract$confidence_levels <- c(
      absent = -1L, unknown = 0L, low = 1L, medium = 2L, high = 3L
    )
    answer$completion_contract$redundant_path_search <- list(
      maximum_paths = options$redundancies,
      cost_increase = options$cost_increase,
      absent_penalty_factor = options$penalty_factor,
      support_threshold = options$support,
      target_flux = options$target_flux
    )
    answer$completion_contract$intentional_corrections <-
      options$intentional_corrections
    answer$union_gem_policy <- paste(
      "one corrected Python-CORDA2 reconstruction per cell type and medium;",
      "five directional confidence levels, redundant-path cost increases,",
      "absent support counting, independent medium testing and final free",
      "reaction completion; only included HC core directions are scored"
    )
    answer$method <- paste(
      "microCOMPASS directional LP on cell-type-specific medium models",
      "reconstructed by corrected Python CORDA2"
    )
    rc_export_microcompass(answer, outdir)
    model_dir <- file.path(outdir, "03_models")
    dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(
      answer$completion_contract,
      file.path(model_dir, "model_completion_contract.rds")
    )
    saveRDS(answer, file.path(outdir, "step_layer2.rds"))
  }
  answer
}
