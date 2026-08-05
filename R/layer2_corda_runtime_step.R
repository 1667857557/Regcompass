# CORDA2 preparation and completion helpers used directly by Layer 2.

.rc_layer2_prepare_completion <- function(
    layer1, meta_modules, model_mode, layer2_args) {
  if (!is.list(layer2_args)) {
    stop("`layer2_args` must be a list.", call. = FALSE)
  }
  model_params <- layer2_args$model_params %||% list()
  if (!is.list(model_params)) {
    stop("`layer2_args$model_params` must be a list.", call. = FALSE)
  }
  is_full_gem <- identical(model_mode, "full_gem")

  if (is_full_gem) {
    .rc_validate_full_gem_model_params(model_params)
    corda_options <- list(
      model_completion = "none",
      requested_model_completion = as.character(
        model_params$model_completion %||% "none"
      ),
      algorithm = "medium_flux_consistency_pruned_full_gem"
    )
    is_corda2 <- FALSE
  } else {
    corda_options <- .rc_layer2_corda_options(model_params)
    is_corda2 <- .rc_is_corda2_options(corda_options)
  }

  if (isTRUE(is_corda2) && !identical(model_mode, "meta_module_gem")) {
    stop(
      "`model_completion = \"corda2\"` is available only with ",
      "`model_mode = \"meta_module_gem\"`.",
      call. = FALSE
    )
  }

  extracted <- c(
    "model_completion", "corda2_args",
    "corda_medium_confidence_threshold",
    "corda_negative_confidence_threshold",
    "corda_regulatory_weight",
    "corda_include_evidence_outside_modules",
    "corda_max_medium_confidence_reactions"
  )
  clean_params <- model_params[setdiff(names(model_params), extracted)]
  layer2_args$model_params <- clean_params

  previous <- as.list(.rc_layer2_completion_context)
  .rc_layer2_completion_context$active <- TRUE
  .rc_layer2_completion_context$model_completion <- if (is_full_gem) {
    "none"
  } else if (isTRUE(is_corda2)) {
    "corda2"
  } else {
    "fastcore"
  }
  .rc_layer2_completion_context$corda_options <- corda_options
  .rc_layer2_completion_context$solver <-
    as.character(layer2_args$solver %||% "highs")
  .rc_layer2_completion_context$completion_time_limit <-
    as.numeric(model_params$completion_time_limit %||% 300)
  .rc_layer2_completion_context$flux_threshold <-
    as.numeric(layer2_args$flux_threshold %||% 1e-8)
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
    is_full_gem = is_full_gem,
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
  if (identical(answer$model_mode, "full_gem")) {
    summary <- answer$model_cache_summary
    contract_summary <- if (is.data.frame(summary)) {
      keep <- intersect(c(
        "medium_scenario", "condition", "n_input_reactions",
        "n_reactions", "n_flux_inconsistent_reactions",
        "flux_consistency_epsilon", "solver",
        "completion_time_limit", "build_strategy"
      ), colnames(summary))
      summary[, keep, drop = FALSE]
    } else {
      data.frame()
    }
    answer$params$model_completion <- "none"
    answer$params$structural_completion <- "medium_flux_consistency"
    answer$params$structural_completion_algorithm <-
      "medium_flux_consistency_pruned_full_gem"
    answer$params$fastcore_executed <- FALSE
    answer$params$corda2_executed <- FALSE
    answer$completion_contract <- list(
      model_completion = "none",
      default_unchanged = FALSE,
      algorithm = "medium_flux_consistency_pruned_full_gem",
      context_specific_reconstruction = FALSE,
      fastcore_executed = FALSE,
      corda2_executed = FALSE,
      medium_applied = TRUE,
      flux_consistency_pruning = TRUE,
      flux_consistency_algorithm = "FASTCC_flux_consistency_only",
      reaction_evidence_used_for_structure = FALSE,
      model_summary = contract_summary
    )
    answer$method <-
      "microCOMPASS shared medium-pruned full-GEM directional LP"
    return(answer)
  }

  answer$params$model_completion <- if (isTRUE(is_corda2)) {
    "corda2"
  } else {
    "fastcore"
  }
  if (!isTRUE(is_corda2)) {
    answer$params$structural_completion <- "fastcore"
    answer$params$structural_completion_algorithm <-
      "add_only_compact_FASTCORE"
    answer$params$fastcore_executed <- TRUE
    answer$params$corda2_executed <- FALSE
    answer$completion_contract <- list(
      model_completion = "fastcore",
      default_unchanged = TRUE,
      algorithm = "add_only_compact_FASTCORE",
      fastcore_executed = TRUE,
      corda2_executed = FALSE
    )
    return(answer)
  }

  original_args <- list(
    MCxNCthresh = corda_options$MCxNCthresh,
    constraint = corda_options$constraint,
    constrainby = corda_options$constrainby,
    om = corda_options$om,
    ci = corda_options$ci
  )
  answer$completion_contract <- list(
    model_completion = "corda2",
    default_unchanged = FALSE,
    algorithm = corda_options$algorithm,
    source_fidelity = "original_MATLAB_CORDA2",
    reference = list(
      repository = corda_options$reference_repository,
      file = corda_options$reference_file
    ),
    adjustable_args = original_args,
    confidence_levels = c(HC = 3L, MC = 2L, NC = 1L, OT = 0L),
    fixed_internal = list(
      fluxThreshold = corda_options$flux_threshold,
      baselineCost = corda_options$baseline_cost,
      outputBound = corda_options$output_bound
    ),
    solver_configuration = list(
      solver = solver,
      threads = 1L
    ),
    stage_update_policy = "original_matlab_directional_order",
    target_parallelism = FALSE,
    fastcore_executed = FALSE,
    corda2_executed = TRUE,
    options = corda_options
  )
  answer$params$structural_completion <- "corda2"
  answer$params$fastcore_executed <- FALSE
  answer$params$corda2_executed <- TRUE
  answer$params$structural_completion_algorithm <- corda_options$algorithm
  answer$params$corda2_args <- original_args
  answer$params$corda2_MCxNCthresh <- corda_options$MCxNCthresh
  answer$params$corda2_constraint <- corda_options$constraint
  answer$params$corda2_constrainby <- corda_options$constrainby
  answer$params$corda2_om <- corda_options$om
  answer$params$corda2_ci <- corda_options$ci
  answer$params$corda2_inner_target_parallelism <- FALSE
  answer$union_gem_policy <- paste(
    "one original-CORDA2 reconstruction per cell type and medium;",
    "targets are processed serially and independent models may run in parallel"
  )
  answer$method <- paste(
    "microCOMPASS directional LP on cell-type-specific medium models",
    "reconstructed with original MATLAB CORDA2 semantics"
  )
  answer
}
