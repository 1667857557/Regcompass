# Exact CORDA2 cell-type by medium cache construction helper.

.rc_build_celltype_medium_corda_cache <- function(
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
      !identical(context$model_completion, "corda2") ||
      !.rc_is_corda2_options(context$corda_options)) {
    stop("CORDA2 cache construction was requested without active options.",
         call. = FALSE)
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
    stop("CORDA2 reaction evidence is unavailable.", call. = FALSE)
  }
  medium_scenarios <- .rc_validate_shared_medium(medium_scenarios)
  scenarios <- unique(as.character(medium_scenarios$medium_scenario_id))
  scenarios <- scenarios[!is.na(scenarios) & nzchar(scenarios)]
  if (!length(scenarios)) {
    stop("No medium scenarios are available for CORDA2 construction.",
         call. = FALSE)
  }
  cache_dir <- file.path(cache_dir, "corda2")
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
    performance <- model$corda_reconstruction$solver_performance %||% list()
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
      n_corda2_lp_solves = as.integer(performance$n_solves %||% NA_integer_),
      n_corda2_objective_coeff_updates = as.integer(
        performance$n_objective_coeff_updates %||% NA_integer_
      ),
      n_corda2_bound_index_updates = as.integer(
        performance$n_bound_index_updates %||% NA_integer_
      ),
      n_corda2_full_vector_values_avoided = as.numeric(
        performance$n_full_vector_numeric_values_avoided %||% NA_real_
      ),
      corda2_transmitted_fraction_of_full = as.numeric(
        performance$transmitted_fraction_of_full %||% NA_real_
      ),
      corda2_solver_release_policy = as.character(
        performance$release_policy %||% NA_character_
      ),
      solver_runtime = paste(execution_types, collapse = ";"),
      target_status = model$target_status,
      build_strategy = "celltype_medium_python_corda2_exact",
      completion_stage = "python_CORDA2_exact_after_confidence_mapping",
      completion_method = "corda2",
      completion_time_limit = Inf,
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
          build_strategy = "celltype_medium_python_corda2_exact"
        )
      }
    }

    rm(model, membership, core, evidence, medium, build, performance)
    list(cache = cache, summary = summary)
  }

  task_bpparam <- .rc_corda_tune_task_bpparam(
    .rc_layer2_task_bpparam(), length(tasks)
  )
  pool_workers <- .rc_corda_pool_workers(task_bpparam)
  outer_parallel <- .rc_corda_should_outer_parallel(
    length(tasks), pool_workers
  )
  active_workers <- if (outer_parallel) {
    min(length(tasks), pool_workers)
  } else {
    1L
  }

  if (outer_parallel) {
    parts <- rc_parallel_lapply(
      tasks,
      function(task) run_one(
        task[1, , drop = FALSE], suppress_nested = TRUE
      ),
      BPPARAM = task_bpparam
    )
    dispatch <-
      "dynamic_cell_type_x_medium_outer_parallel_python_instances"
  } else {
    parts <- lapply(tasks, function(task) {
      run_one(task[1, , drop = FALSE], suppress_nested = FALSE)
    })
    dispatch <- "serial_within_each_python_corda2_instance"
  }

  cache <- list()
  summaries <- vector("list", length(parts))
  for (i in seq_along(parts)) {
    part <- parts[[i]]
    summaries[[i]] <- part$summary
    if (length(part$cache)) {
      duplicated <- intersect(names(cache), names(part$cache))
      if (length(duplicated)) {
        stop("Parallel CORDA2 tasks produced duplicate cache keys.",
             call. = FALSE)
      }
      cache[names(part$cache)] <- part$cache
    }
  }
  attr(cache, "summary") <- .rc_bind_frames_fill(summaries)
  attr(cache, "celltype_col") <- celltype_col
  attr(cache, "structural_scope") <- "cell_type_x_medium"
  attr(cache, "completion_method") <- "corda2"
  attr(cache, "structural_parallel_task") <- dispatch
  attr(cache, "structural_parallel_workers_requested") <- pool_workers
  attr(cache, "structural_parallel_workers") <- active_workers
  attr(cache, "structural_parallel_tasks") <- length(tasks)
  attr(cache, "structural_dynamic_task_scheduling") <- outer_parallel
  attr(cache, "corda2_inner_target_parallelism") <- FALSE
  attr(cache, "fastcore_parallel_task") <- "not_applicable_to_corda2"
  rm(parts, summaries)
  invisible(gc(verbose = FALSE, full = FALSE))
  cache
}
