# Runtime dispatch for optional original CORDA Layer 2 completion.

.rc_layer2_completion_context <- new.env(parent = emptyenv())
.rc_layer2_completion_context$active <- FALSE
.rc_layer2_completion_context$model_completion <- "fastcore"
.rc_layer2_completion_context$reaction_evidence <- NULL
.rc_layer2_completion_context$corda_options <- NULL

.rc_regcompass_step_layer2_completion_base <- rc_regcompass_step_layer2
.rc_build_celltype_medium_union_gem_cache_fastcore <-
  .rc_build_celltype_medium_union_gem_cache

# Correct persistent-solver fallback while retaining basis reuse on success.
.rc_corda_engine_solve <- function(
    engine, objective, lower, upper) {
  engine$n_solves <- engine$n_solves + 1L
  if (!identical(engine$type, "highs_persistent_cpp")) {
    return(list(
      answer = .rc_corda_one_shot_solve(engine, objective, lower, upper),
      engine = engine
    ))
  }
  persistent <- tryCatch({
    index <- seq_along(objective) - 1L
    .rc_corda_highs_call(
      "hi_solver_set_objective", engine$pointer,
      index = index, coeff = as.numeric(objective)
    )
    .rc_corda_highs_call(
      "hi_solver_set_variable_bounds", engine$pointer,
      index = index,
      lower = as.numeric(lower),
      upper = as.numeric(upper)
    )
    .rc_corda_highs_call("hi_solver_run", engine$pointer)
    status_message <- .rc_corda_highs_call(
      "hi_solver_status_message", engine$pointer
    )
    status <- .rc_lp_status(status_message)
    solution <- if (identical(status, "optimal")) {
      as.numeric(.rc_corda_highs_call(
        "hi_solver_get_solution", engine$pointer
      )$col_value)
    } else {
      numeric()
    }
    info <- tryCatch(
      .rc_corda_highs_call("hi_solver_info", engine$pointer),
      error = function(e) list()
    )
    list(
      status = status,
      solution = solution,
      objective = as.numeric(
        info$objective_function_value %||% NA_real_
      ),
      backend = "highs_persistent_cpp_basis_reuse",
      solver_message = as.character(status_message)
    )
  }, error = function(e) e)
  if (!inherits(persistent, "error")) {
    return(list(answer = persistent, engine = engine))
  }
  failure_message <- conditionMessage(persistent)
  engine$n_fallback <- engine$n_fallback + 1L
  fallback <- .rc_corda_one_shot_solve(engine, objective, lower, upper)
  fallback$backend <- "highs_persistent_failed_one_shot_fallback"
  fallback$solver_message <- paste(
    failure_message, fallback$solver_message
  )
  list(answer = fallback, engine = engine)
}

.rc_corda_pool_workers <- function(BPPARAM) {
  if (identical(BPPARAM, FALSE)) return(1L)
  if (!is.null(BPPARAM) &&
      requireNamespace("BiocParallel", quietly = TRUE) &&
      methods::is(BPPARAM, "BiocParallelParam")) {
    return(max(1L, BiocParallel::bpnworkers(BPPARAM)))
  }
  config <- rc_parallel_config(workers = NULL, backend = "auto")
  max(1L, config$workers)
}

.rc_build_celltype_medium_union_gem_cache <- function(
    gem, reaction_membership, core_reactions,
    target_reactions = NULL, medium_scenarios = NULL,
    celltype_col = "cell_type",
    cache_dir = tempfile("RegCompassR_celltype_medium_union_gem_"),
    target_direction = c("both", "forward", "reverse"),
    solver = "highs", time_limit = 300,
    fastcore_epsilon = 1e-4,
    max_support_reactions = 2000,
    strict = TRUE) {
  context <- .rc_layer2_completion_context
  if (!isTRUE(context$active) ||
      !identical(context$model_completion, "corda")) {
    return(.rc_build_celltype_medium_union_gem_cache_fastcore(
      gem = gem,
      reaction_membership = reaction_membership,
      core_reactions = core_reactions,
      target_reactions = target_reactions,
      medium_scenarios = medium_scenarios,
      celltype_col = celltype_col,
      cache_dir = cache_dir,
      target_direction = target_direction,
      solver = solver,
      time_limit = time_limit,
      fastcore_epsilon = fastcore_epsilon,
      max_support_reactions = max_support_reactions,
      strict = strict
    ))
  }
  target_direction <- match.arg(target_direction)
  scoped <- .rc_validate_celltype_reaction_scope(
    reaction_membership, core_reactions, celltype_col
  )
  reaction_membership <- scoped$reaction_membership
  core_reactions <- scoped$core_reactions
  evidence_all <- context$reaction_evidence
  corda_options <- context$corda_options
  if (!is.data.frame(evidence_all) || !nrow(evidence_all)) {
    stop("CORDA reaction evidence is unavailable.", call. = FALSE)
  }
  medium_scenarios <- .rc_validate_shared_medium(medium_scenarios)
  scenarios <- unique(as.character(medium_scenarios$medium_scenario_id))
  scenarios <- scenarios[!is.na(scenarios) & nzchar(scenarios)]
  if (!length(scenarios)) {
    stop("No medium scenarios are available for CORDA construction.",
         call. = FALSE)
  }
  cache_dir <- file.path(cache_dir, "corda")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  task_grid <- expand.grid(
    cell_type = scoped$cell_types,
    medium_scenario = scenarios,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  tasks <- split(task_grid, seq_len(nrow(task_grid)))
  run_one <- function(task, suppress_nested = FALSE) {
    previous_nested <- .rc_layer2_parallel_context$nested_serial
    if (isTRUE(suppress_nested)) {
      .rc_layer2_parallel_context$nested_serial <- TRUE
    }
    on.exit({
      .rc_layer2_parallel_context$nested_serial <- previous_nested
      invisible(gc(verbose = FALSE, full = TRUE))
    }, add = TRUE)
    cell_type <- as.character(task$cell_type[[1L]])
    scenario <- as.character(task$medium_scenario[[1L]])
    membership <- reaction_membership[
      reaction_membership[[celltype_col]] == cell_type, , drop = FALSE
    ]
    core <- core_reactions[
      core_reactions[[celltype_col]] == cell_type, , drop = FALSE
    ]
    selected_targets <- .rc_celltype_target_reactions(
      target_reactions, cell_type, celltype_col
    )
    if (!is.null(selected_targets)) {
      core <- core[core$reaction_id %in% selected_targets, , drop = FALSE]
    }
    if (!nrow(core)) {
      stop("No core targets remain for cell type `", cell_type, "`.",
           call. = FALSE)
    }
    medium <- medium_scenarios[
      as.character(medium_scenarios$medium_scenario_id) == scenario,
      , drop = FALSE
    ]
    if (!nrow(medium) ||
        (".no_constraints" %in% colnames(medium) &&
         all(medium$.no_constraints))) {
      medium <- NULL
    }
    evidence <- evidence_all[
      as.character(evidence_all$cell_type) == cell_type, , drop = FALSE
    ]
    model <- .rc_complete_celltype_medium_corda_gem(
      gem = gem,
      reaction_membership = membership,
      core_reactions = core,
      cell_type = cell_type,
      reaction_evidence = evidence,
      corda_options = corda_options,
      medium_table = medium,
      target_direction = target_direction,
      solver = solver,
      time_limit = time_limit,
      fastcore_epsilon = fastcore_epsilon,
      max_support_reactions = max_support_reactions,
      strict = strict
    )
    model$shared_across_conditions <- TRUE
    model$shared_across_cell_types <- FALSE
    model$is_union_gem <- TRUE
    model$union_gem_cell_type <- cell_type
    model$union_gem_medium_scenario <- scenario
    model$union_gem_scope <-
      "one_cell_type_one_medium_shared_across_conditions_and_matching_metacells"
    file <- file.path(
      cache_dir,
      paste0(
        "union_gem__celltype_", .rc_safe_cache_token(cell_type),
        "__medium_", .rc_safe_cache_token(scenario), ".rds"
      )
    )
    .rc_atomic_save_rds(model, file)
    checksum <- unname(tools::md5sum(file))
    build <- model$build_params
    execution_types <- unique(unlist(lapply(
      model$corda_execution,
      function(value) value$solver_runtime %||% character()
    ), use.names = FALSE))
    summary <- data.frame(
      cell_type = cell_type,
      medium_scenario = scenario,
      file = file,
      file_checksum = checksum,
      n_reactions = ncol(model$S),
      n_metabolites = nrow(model$S),
      n_celltype_biological_reactions =
        build$n_celltype_biological_reactions,
      n_celltype_fastcore_support_reactions =
        build$n_celltype_fastcore_support_reactions,
      n_high_confidence_reactions = build$n_high_confidence_reactions,
      n_module_medium_confidence_reactions =
        build$n_module_medium_confidence_reactions,
      n_evidence_medium_confidence_reactions =
        build$n_evidence_medium_confidence_reactions,
      n_negative_confidence_reactions =
        build$n_negative_confidence_reactions,
      n_other_reactions = build$n_other_reactions,
      n_corda_included_reactions = build$n_corda_included_reactions,
      n_corda_included_initial_MC = build$n_corda_included_initial_MC,
      n_corda_included_initial_NC = build$n_corda_included_initial_NC,
      n_corda_included_initial_OT = build$n_corda_included_initial_OT,
      n_stage1_associated = build$n_stage1_associated,
      n_stage2_promoted_nc = build$n_stage2_promoted_nc,
      n_stage2_promoted_mc = build$n_stage2_promoted_mc,
      n_stage3_associated_ot = build$n_stage3_associated_ot,
      solver_runtime = paste(execution_types, collapse = ";"),
      target_status = model$target_status,
      build_strategy = "celltype_medium_original_corda",
      completion_stage =
        "original_CORDA_three_stage_after_confidence_mapping",
      completion_time_limit = time_limit,
      stringsAsFactors = FALSE
    )
    cache <- list()
    if (nrow(model$target_directions)) {
      for (i in seq_len(nrow(model$target_directions))) {
        reaction <- as.character(model$target_directions$reaction_id[[i]])
        direction <- as.character(
          model$target_directions$target_direction[[i]]
        )
        key <- paste0(
          "celltype=", utils::URLencode(cell_type, reserved = TRUE),
          "::reaction=", utils::URLencode(reaction, reserved = TRUE),
          "::direction=", direction,
          "::medium=", utils::URLencode(scenario, reserved = TRUE)
        )
        cache[[key]] <- list(
          module_id = "CELLTYPE_MEDIUM_UNION_GEM",
          cell_type = cell_type,
          reaction_id = reaction,
          target_direction = direction,
          medium_scenario = scenario,
          condition = "all",
          file = file,
          file_checksum = checksum,
          build_strategy = "celltype_medium_original_corda"
        )
      }
    }
    rm(model)
    list(cache = cache, summary = summary)
  }

  task_bpparam <- .rc_layer2_task_bpparam()
  pool_workers <- .rc_corda_pool_workers(task_bpparam)
  outer_parallel <- length(tasks) > 1L &&
    length(tasks) >= pool_workers && pool_workers > 1L
  if (outer_parallel) {
    parts <- rc_parallel_lapply(
      tasks,
      function(task) run_one(
        task[1, , drop = FALSE], suppress_nested = TRUE
      ),
      BPPARAM = task_bpparam
    )
    dispatch <- "cell_type_x_medium_outer_parallel"
  } else {
    parts <- lapply(tasks, function(task) {
      run_one(task[1, , drop = FALSE], suppress_nested = FALSE)
    })
    dispatch <- if (pool_workers > 1L) {
      "models_serial_target_direction_x_replicate_inner_parallel"
    } else {
      "serial"
    }
  }
  cache <- list()
  summaries <- vector("list", length(parts))
  for (i in seq_along(parts)) {
    part <- parts[[i]]
    summaries[[i]] <- part$summary
    if (length(part$cache)) {
      duplicated <- intersect(names(cache), names(part$cache))
      if (length(duplicated)) {
        stop("Parallel CORDA tasks produced duplicate cache keys.",
             call. = FALSE)
      }
      cache[names(part$cache)] <- part$cache
    }
  }
  attr(cache, "summary") <- .rc_bind_frames_fill(summaries)
  attr(cache, "celltype_col") <- celltype_col
  attr(cache, "structural_scope") <- "cell_type_x_medium"
  attr(cache, "completion_method") <- "corda"
  attr(cache, "structural_parallel_task") <- dispatch
  attr(cache, "structural_parallel_workers") <- pool_workers
  attr(cache, "fastcore_parallel_task") <- "not_applicable_to_corda"
  cache
}

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
  if (identical(corda_options$model_completion, "corda") &&
      !identical(model_mode, "meta_module_gem")) {
    stop(
      "`model_completion = \"corda\"` is available only with ",
      "`model_mode = \"meta_module_gem\"`.",
      call. = FALSE
    )
  }
  extracted <- c(
    "model_completion",
    "corda_gamma", "corda_kappa", "corda_epsilon",
    "corda_n", "corda_p", "corda_seed", "corda_flux_tolerance",
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
    corda_options$model_completion
  .rc_layer2_completion_context$corda_options <- corda_options
  .rc_layer2_completion_context$reaction_evidence <- if (
    identical(corda_options$model_completion, "corda")
  ) {
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
  answer$params$model_completion <- corda_options$model_completion
  answer$completion_contract <- list(
    model_completion = corda_options$model_completion,
    default_unchanged = identical(corda_options$model_completion, "fastcore"),
    algorithm = if (identical(corda_options$model_completion, "corda")) {
      "Schultz_Qutub_CORDA_2016_three_stage_dependency_assessment"
    } else {
      "add_only_compact_FASTCORE"
    },
    confidence_mapping = if (
      identical(corda_options$model_completion, "corda")
    ) {
      paste(
        "HC=merged core; MC=remaining module plus optional high-evidence",
        "outside-module reactions; NC=low evidence; OT=remaining"
      )
    } else {
      NULL
    },
    stage_update_policy = if (
      identical(corda_options$model_completion, "corda")
    ) "barrier_then_union_order_independent" else NULL,
    native_solver_acceleration = if (
      identical(corda_options$model_completion, "corda")
    ) {
      "persistent HiGHS C++ solver with objective/bound updates and basis reuse"
    } else {
      NULL
    },
    corda = if (
      identical(corda_options$model_completion, "corda")
    ) corda_options else NULL
  )
  if (identical(corda_options$model_completion, "corda")) {
    answer$union_gem_policy <- paste(
      "one original-CORDA reconstruction per cell type and medium;",
      "three dependency-assessment stages use barrier updates; only retained",
      "parent-feasible HC core directions enter downstream scoring"
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
