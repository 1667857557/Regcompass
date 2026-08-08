# CORDA2 preparation and completion helpers used directly by Layer 2.

.rc_layer2_completion_time_limit <- function(model_params, is_corda2) {
  if (isTRUE(is_corda2)) {
    if ("completion_time_limit" %in% names(model_params)) {
      stop(
        "`layer2_args$model_params$completion_time_limit` is not supported ",
        "for CORDA2. CORDA2 reconstruction runs without a time limit; ",
        "remove this parameter.",
        call. = FALSE
      )
    }
    return(Inf)
  }

  value <- suppressWarnings(as.numeric(
    model_params$completion_time_limit %||% 300
  ))
  if (length(value) != 1L || is.na(value) || value <= 0) {
    stop(
      "`layer2_args$model_params$completion_time_limit` must be a positive ",
      "number or Inf.",
      call. = FALSE
    )
  }
  value
}

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
      algorithm = "compass_medium_constrained_full_gem"
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

  completion_time_limit <- .rc_layer2_completion_time_limit(
    model_params, is_corda2
  )

  extracted <- c(
    "model_completion", "corda2_args",
    "corda_medium_confidence_threshold",
    "corda_negative_confidence_threshold",
    "corda_regulatory_weight",
    "corda_include_evidence_outside_modules",
    "corda_max_medium_confidence_reactions"
  )
  clean_params <- model_params[setdiff(names(model_params), extracted)]
  if (isTRUE(is_corda2)) {
    clean_params$completion_time_limit <- Inf
  }
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
    completion_time_limit
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
        "n_reactions", "n_medium_removed_reactions",
        "n_medium_bound_changes", "medium_applied",
        "medium_handling", "build_strategy"
      ), colnames(summary))
      summary[, keep, drop = FALSE]
    } else {
      data.frame()
    }
    answer$params$model_completion <- "none"
    answer$params$structural_completion <- "none"
    answer$params$structural_completion_algorithm <-
      "compass_medium_bounds_only"
    answer$params$medium_handling <-
      "exchange_bounds_only_no_reaction_deletion"
    answer$params$medium_direct_reaction_deletion <- FALSE
    answer$params$fastcore_executed <- FALSE
    answer$params$corda2_executed <- FALSE
    answer$completion_contract <- list(
      model_completion = "none",
      default_unchanged = FALSE,
      algorithm = "compass_medium_constrained_full_gem",
      context_specific_reconstruction = FALSE,
      fastcore_executed = FALSE,
      corda2_executed = FALSE,
      medium_applied = if (
        is.data.frame(summary) && "medium_applied" %in% colnames(summary)
      ) all(summary$medium_applied %in% TRUE) else NA,
      medium_handling = "exchange_bounds_only_no_reaction_deletion",
      medium_direct_reaction_deletion = FALSE,
      flux_consistency_pruning = FALSE,
      flux_consistency_algorithm = "none",
      target_feasibility = paste(
        "retain the complete medium-constrained GEM; compute directional vmax;",
        "skip the penalty LP when vmax is below the scoring threshold"
      ),
      reaction_evidence_used_for_structure = FALSE,
      model_summary = contract_summary
    )
    answer$method <-
      "microCOMPASS shared medium-constrained full-GEM directional LP"
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
    answer$params$medium_handling <-
      "exchange_bounds_only_then_fastcc_fastcore"
    answer$params$medium_direct_reaction_deletion <- FALSE
    answer$params$fastcore_executed <- TRUE
    answer$params$corda2_executed <- FALSE
    answer$completion_contract <- list(
      model_completion = "fastcore",
      default_unchanged = FALSE,
      algorithm = "add_only_compact_FASTCORE",
      medium_handling = "exchange_bounds_only_then_fastcc_fastcore",
      medium_direct_reaction_deletion = FALSE,
      fastcc_role = paste(
        "FASTCC is part of FASTCORE parent consistency analysis;",
        "the medium table itself only changes exchange bounds"
      ),
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
  stage_parallel <- isTRUE(attr(
    answer$shared_model_cache, "corda2_inner_target_parallelism"
  ))
  execution_policy <- if (stage_parallel) {
    "stage_barrier_parallel_targets_deterministic_ordered_reduce"
  } else {
    "serial_original_persistent_engine"
  }
  target_parallelism <- if (stage_parallel) {
    "within_each_corda2_stage"
  } else {
    FALSE
  }
  worker_lifecycle <- if (stage_parallel) {
    "fresh_pool_each_stage_stop_full_gc_before_next_stage"
  } else {
    "single_persistent_engine_for_complete_reconstruction"
  }
  answer$completion_contract <- list(
    model_completion = "corda2",
    default_unchanged = TRUE,
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
      threads = 1L,
      completion_time_limit = Inf
    ),
    stage_update_policy = "original_matlab_directional_order",
    parallel_execution_policy = execution_policy,
    target_parallelism = target_parallelism,
    stage_barrier = stage_parallel,
    stage_worker_lifecycle = worker_lifecycle,
    medium_handling = "exchange_bounds_only_then_corda2",
    medium_direct_reaction_deletion = FALSE,
    parent_prepruning = "none",
    fastcore_executed = FALSE,
    corda2_executed = TRUE,
    options = corda_options
  )
  answer$params$structural_completion <- "corda2"
  answer$params$medium_handling <- "exchange_bounds_only_then_corda2"
  answer$params$medium_direct_reaction_deletion <- FALSE
  answer$params$fastcore_executed <- FALSE
  answer$params$corda2_executed <- TRUE
  answer$params$structural_completion_algorithm <- corda_options$algorithm
  answer$params$corda2_args <- original_args
  answer$params$corda2_MCxNCthresh <- corda_options$MCxNCthresh
  answer$params$corda2_constraint <- corda_options$constraint
  answer$params$corda2_constrainby <- corda_options$constrainby
  answer$params$corda2_om <- corda_options$om
  answer$params$corda2_ci <- corda_options$ci
  answer$params$corda2_completion_time_limit <- Inf
  answer$params$corda2_inner_target_parallelism <- stage_parallel
  answer$params$corda2_stage_barrier_parallelism <- stage_parallel
  answer$union_gem_policy <- if (stage_parallel) {
    paste(
      "one original-CORDA2 reconstruction per cell type and medium;",
      "Step 1, Step 2.1, Step 2.2 and Step 3 remain strict barriers;",
      "directional targets inside each step use the full Layer-2 worker budget"
    )
  } else {
    paste(
      "one original-CORDA2 reconstruction per cell type and medium;",
      "the preserved original persistent-engine serial target order is used"
    )
  }
  answer$method <- paste(
    "microCOMPASS directional LP on cell-type-specific medium models",
    "reconstructed with original MATLAB CORDA2 semantics"
  )
  answer
}
