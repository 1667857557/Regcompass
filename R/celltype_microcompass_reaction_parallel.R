# Reaction-granular directional LP scoring on cell-type-specific models.

.rc_run_celltype_microcompass_engine_reaction_core <- function(
    layer1, gem, target_reactions = NULL,
    medium_table = NULL, medium_scenarios = NULL,
    reaction_membership, core_reactions,
    unit = c("metacell", "sample_celltype"),
    condition_col = "condition", sample_col = NULL,
    celltype_col = "cell_type", model_params = list(),
    omega = 0.95,
    target_direction = c("both", "forward", "reverse"),
    parallel = TRUE,
    solver = c("highs", "gurobi", "glpk"),
    flux_threshold = 1e-8,
    BPPARAM = NULL,
    model_cache_override = NULL) {
  unit <- match.arg(unit)
  solver <- match.arg(solver)
  target_direction <- match.arg(target_direction)
  medium_scenarios <- .rc_validate_shared_medium(
    medium_scenarios %||% medium_table
  )
  matrices <- rc_layer2_unit_matrices(
    layer1,
    if (identical(unit, "metacell")) "metacell" else "sample_celltype",
    sample_col, celltype_col, condition_col
  )
  unit_meta <- matrices$unit_meta
  if (!is.data.frame(unit_meta) || !celltype_col %in% colnames(unit_meta)) {
    stop("Layer 1 unit metadata lack the cell-type column.", call. = FALSE)
  }
  unit_id <- if ("unit_id" %in% colnames(unit_meta)) {
    as.character(unit_meta$unit_id)
  } else if ("pool_id" %in% colnames(unit_meta)) {
    as.character(unit_meta$pool_id)
  } else {
    stop("Layer 1 unit metadata lack unit_id/pool_id.", call. = FALSE)
  }
  if (!setequal(unit_id, colnames(matrices$reaction_expression))) {
    stop("Layer 1 unit metadata and expression columns differ.", call. = FALSE)
  }
  unit_meta <- unit_meta[
    match(colnames(matrices$reaction_expression), unit_id), , drop = FALSE
  ]
  unit_celltype <- trimws(as.character(unit_meta[[celltype_col]]))
  if (anyNA(unit_celltype) || any(!nzchar(unit_celltype))) {
    stop("Layer 1 unit cell types must be complete.", call. = FALSE)
  }
  names(unit_celltype) <- colnames(matrices$reaction_expression)
  gem <- rc_annotate_reaction_roles(gem)

  if (!is.null(model_cache_override)) {
    if (!is.list(model_cache_override) || !length(model_cache_override) ||
        is.null(attr(model_cache_override, "summary"))) {
      stop("`model_cache_override` is not an audited cell-type model cache.",
           call. = FALSE)
    }
    model_cache <- model_cache_override
  } else {
    model_cache <- .rc_build_celltype_medium_union_gem_cache(
      gem = gem,
      reaction_membership = reaction_membership,
      core_reactions = core_reactions,
      target_reactions = target_reactions,
      medium_scenarios = medium_scenarios,
      celltype_col = celltype_col,
      cache_dir = model_params$cache_dir %||%
        tempfile("RegCompassR_celltype_union_gem_cache_"),
      target_direction = target_direction,
      solver = solver,
      time_limit = model_params$completion_time_limit %||% 300,
      fastcore_epsilon = model_params$fastcore_epsilon %||% 1e-4,
      max_support_reactions = model_params$max_support_reactions %||% 2000,
      strict = model_params$strict %||% TRUE
    )
  }
  if (!length(model_cache)) {
    stop("No cell-type union-GEM targets were available.", call. = FALSE)
  }
  cache_celltypes <- sort(unique(vapply(
    model_cache, function(entry) as.character(entry$cell_type), character(1)
  )))
  unit_celltypes <- sort(unique(unit_celltype))
  if (is.null(model_cache_override)) {
    if (!setequal(cache_celltypes, unit_celltypes)) {
      stop(
        "Stage 5 union GEMs and Layer 1 units cover different cell types.",
        call. = FALSE
      )
    }
  } else if (!all(cache_celltypes %in% unit_celltypes)) {
    stop(
      "A reused union-GEM cache contains cell types absent from Layer 1.",
      call. = FALSE
    )
  }
  summary <- attr(model_cache, "summary")
  requested_media <- unique(as.character(medium_scenarios$medium_scenario_id))
  if (!is.data.frame(summary) ||
      !all(c("cell_type", "medium_scenario") %in% colnames(summary)) ||
      !setequal(unique(as.character(summary$medium_scenario)), requested_media)) {
    stop("Cell-type model cache and requested media differ.", call. = FALSE)
  }

  row_ids <- names(model_cache)
  units <- colnames(matrices$reaction_expression)
  model_keys <- vapply(model_cache, function(entry) entry$file, character(1))
  names(model_keys) <- row_ids
  unique_model_keys <- unique(unname(model_keys))
  representative_rows <- vapply(
    unique_model_keys,
    function(key) names(model_keys)[match(key, model_keys)],
    character(1)
  )
  names(representative_rows) <- unique_model_keys
  all_reactions <- unique(unlist(lapply(representative_rows, function(row_id) {
    entry <- model_cache[[row_id]]
    colnames(.rc_read_celltype_union_gem(
      entry$file, entry$cell_type, entry$medium_scenario,
      entry$file_checksum
    )$S)
  }), use.names = FALSE))
  penalties <- rc_compute_multiome_penalty(
    rc_align_reaction_expression(
      matrices$reaction_expression, all_reactions, NA_real_
    ),
    reaction_roles = gem$reaction_roles
  )
  vmax_cache <- .rc_build_microcompass_vmax_cache(
    model_cache = model_cache,
    mode = "meta_module_gem",
    model_keys = model_keys,
    solver = solver,
    flux_threshold = flux_threshold,
    parallel = parallel,
    BPPARAM = BPPARAM
  )
  vmax_cache_diagnostics <- do.call(rbind, lapply(row_ids, function(row_id) {
    entry <- model_cache[[row_id]]
    value <- vmax_cache[[row_id]]
    data.frame(
      row_id = row_id,
      cell_type = entry$cell_type,
      reaction_id = entry$reaction_id,
      target_direction = entry$target_direction,
      medium_scenario = entry$medium_scenario,
      vmax = as.numeric(value$vmax),
      feasible = isTRUE(value$feasible),
      status = as.character(value$status),
      computation_scope =
        "directional_target_batch_within_shared_model",
      stringsAsFactors = FALSE
    )
  }))
  rownames(vmax_cache_diagnostics) <- NULL

  penalty <- vmax <- matrix(
    NA_real_, length(row_ids), length(units),
    dimnames = list(row_ids, units)
  )
  feasible <- evaluated <- matrix(
    FALSE, length(row_ids), length(units),
    dimnames = list(row_ids, units)
  )

  checkpoint_root <- model_params$cache_dir %||%
    dirname(model_cache[[row_ids[[1L]]]]$file)
  checkpoint_dir <- tempfile(
    pattern = "step2_reaction_tasks_",
    tmpdir = checkpoint_root
  )
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(checkpoint_dir, recursive = TRUE, force = TRUE), add = TRUE)
  tasks <- stats::setNames(as.list(row_ids), row_ids)

  run_one_reaction <- function(row_id) {
    on.exit(invisible(gc(verbose = FALSE, full = TRUE)), add = TRUE)
    entry <- model_cache[[row_id]]
    eligible <- names(unit_celltype)[unit_celltype == entry$cell_type]
    if (!length(eligible)) {
      stop("No Layer 1 units match cell type `", entry$cell_type, "`.",
           call. = FALSE)
    }
    model <- .rc_read_celltype_union_gem(
      entry$file, entry$cell_type,
      entry$medium_scenario, entry$file_checksum
    )
    target_index <- match(entry$reaction_id, colnames(model$S))
    if (is.na(target_index)) {
      stop("A target reaction is absent from its cell-type union GEM.",
           call. = FALSE)
    }

    prepared <- .rc_compass_step2_prepare(
      S = model$S,
      lb = model$lb,
      ub = model$ub,
      target_reaction = entry$reaction_id,
      vmax_result = vmax_cache[[row_id]],
      target_direction = entry$target_direction,
      omega = omega,
      flux_threshold = flux_threshold
    )
    step2_engine <- if (isTRUE(prepared$runnable)) {
      .rc_compass_step2_new_engine(prepared$template, solver)
    } else {
      NULL
    }
    on.exit(
      .rc_compass_step2_release_engine(step2_engine),
      add = TRUE
    )

    task_penalty <- rep(NA_real_, length(eligible))
    task_vmax <- rep(NA_real_, length(eligible))
    task_feasible <- task_evaluated <- rep(FALSE, length(eligible))
    names(task_penalty) <- names(task_vmax) <-
      names(task_feasible) <- names(task_evaluated) <- eligible
    diagnostics <- vector("list", length(eligible))

    for (i in seq_along(eligible)) {
      one_unit <- eligible[[i]]
      unit_penalty <- penalties$penalty[colnames(model$S), one_unit]
      target_penalty <- unit_penalty[[target_index]]
      evidence_available <- is.finite(unit_penalty)
      solver_penalty <- unit_penalty
      solver_penalty[!evidence_available] <- 0

      fit <- if (isTRUE(prepared$runnable)) {
        solved <- .rc_compass_step2_engine_solve(
          step2_engine, solver_penalty
        )
        step2_engine <- solved$engine
        .rc_compass_step2_result(prepared$template, solved$answer)
      } else {
        .rc_compass_step2_align_penalties(
          prepared$reactions, solver_penalty
        )
        prepared$result
      }

      target_available <- is.finite(target_penalty)
      task_penalty[[one_unit]] <- if (target_available) {
        fit$penalty
      } else {
        NA_real_
      }
      task_vmax[[one_unit]] <- fit$vmax
      task_feasible[[one_unit]] <- isTRUE(fit$feasible)
      task_evaluated[[one_unit]] <-
        isTRUE(fit$feasible) && target_available
      diagnostics[[i]] <- data.frame(
        row_id = row_id,
        unit_id = one_unit,
        module_id = "CELLTYPE_MEDIUM_UNION_GEM",
        cell_type = entry$cell_type,
        reaction_id = entry$reaction_id,
        target_direction = entry$target_direction,
        medium_scenario = entry$medium_scenario,
        condition = "all",
        strict_feasible = isTRUE(fit$feasible),
        solver_status = fit$solver_status,
        solver_backend = fit$solver_backend %||% "unknown",
        step1_status = fit$step1_status,
        step2_status = fit$step2_status,
        target_status = model$target_status %||%
          if (isTRUE(fit$feasible)) "ok" else "structurally_infeasible",
        objective_value = if (target_available) fit$penalty else NA_real_,
        vmax = fit$vmax,
        vmax_reused_from_celltype_cache = TRUE,
        step2_model_reused_across_metacells = TRUE,
        target_expression_available = target_available,
        objective_evidence_fraction = mean(evidence_available),
        unavailable_objective_terms = sum(!evidence_available),
        parallel_task = "directional_reaction_x_matching_metacells",
        stringsAsFactors = FALSE
      )
      rm(unit_penalty, solver_penalty, fit)
    }

    engine_metrics <- .rc_compass_step2_engine_metrics(step2_engine)
    token <- substr(.rc_microcompass_object_checksum(list(
      row_id = row_id,
      file_checksum = entry$file_checksum,
      units = eligible,
      omega = omega,
      solver = solver,
      flux_threshold = flux_threshold
    )), 1L, 24L)
    checkpoint <- file.path(
      checkpoint_dir, paste0("step2__", token, ".rds")
    )
    .rc_atomic_save_rds(list(
      row_id = row_id,
      units = eligible,
      penalty = task_penalty,
      vmax = task_vmax,
      feasible = task_feasible,
      evaluated = task_evaluated,
      diagnostics = .rc_bind_frames_fill(diagnostics),
      engine_metrics = engine_metrics
    ), checkpoint)
    rm(model, diagnostics, task_penalty, task_vmax,
       task_feasible, task_evaluated, prepared)
    checkpoint
  }

  checkpoint_files <- rc_parallel_lapply(
    tasks,
    run_one_reaction,
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )
  diagnostics <- vector("list", length(checkpoint_files))
  step2_engine_metrics <- vector("list", length(checkpoint_files))
  for (i in seq_along(checkpoint_files)) {
    result <- readRDS(checkpoint_files[[i]])
    row_id <- as.character(result$row_id)
    if (!row_id %in% row_ids ||
        !identical(names(result$penalty), as.character(result$units))) {
      stop("A reaction-level Step 2 checkpoint is malformed.",
           call. = FALSE)
    }
    penalty[row_id, result$units] <- result$penalty
    vmax[row_id, result$units] <- result$vmax
    feasible[row_id, result$units] <- result$feasible
    evaluated[row_id, result$units] <- result$evaluated
    diagnostics[[i]] <- result$diagnostics
    metrics <- result$engine_metrics %||% list()
    step2_engine_metrics[[i]] <- data.frame(
      row_id = row_id,
      engine = as.character(metrics$engine %||% "unknown"),
      n_solves = as.integer(metrics$n_solves %||% 0L),
      n_objective_updates = as.integer(
        metrics$n_objective_updates %||% 0L
      ),
      n_fallback = as.integer(metrics$n_fallback %||% 0L),
      stringsAsFactors = FALSE
    )
    rm(result)
    unlink(checkpoint_files[[i]], force = TRUE)
    invisible(gc(verbose = FALSE, full = TRUE))
  }

  score <- rc_compass_score_from_penalty(penalty, feasible)
  directions <- unique(do.call(rbind, lapply(model_cache, function(entry) {
    data.frame(
      cell_type = entry$cell_type,
      reaction_id = entry$reaction_id,
      target_direction = entry$target_direction,
      medium_scenario = entry$medium_scenario,
      stringsAsFactors = FALSE
    )
  })))
  model_diagnostics <- .rc_bind_frames_fill(lapply(
    representative_rows, function(row_id) {
      entry <- model_cache[[row_id]]
      model <- .rc_read_celltype_union_gem(
        entry$file, entry$cell_type, entry$medium_scenario,
        entry$file_checksum
      )
      out <- model$closure_diagnostics %||% data.frame()
      if (nrow(out)) {
        out$cell_type <- entry$cell_type
        out$medium_scenario <- entry$medium_scenario
      }
      out
    }
  ))
  reuse <- table(unit_celltype)
  list(
    score = score,
    penalty = penalty,
    vmax = vmax,
    feasible = feasible,
    evaluated = evaluated,
    target_direction = directions,
    direction_diagnostics = directions,
    medium_scenarios = medium_scenarios,
    model_mode = "meta_module_gem",
    shared_model_cache = model_cache,
    model_cache_summary = summary,
    structural_model_contract = .rc_celltype_model_contract(model_cache),
    model_diagnostics = model_diagnostics,
    vmax_cache_diagnostics = vmax_cache_diagnostics,
    lp_diagnostics = .rc_bind_frames_fill(diagnostics),
    step2_engine_metrics = .rc_bind_frames_fill(step2_engine_metrics),
    penalty_components = penalties$components,
    evidence_policy = penalties$evidence_policy,
    evidence_policy_detail = penalties$evidence_policy_detail,
    unit_meta = unit_meta,
    params = list(
      unit = unit,
      omega = omega,
      target_direction = target_direction,
      shared_gem = TRUE,
      shared_across_conditions = TRUE,
      shared_across_cell_types = FALSE,
      shared_gem_scope =
        "one_union_gem_per_cell_type_per_medium_shared_within_cell_type",
      structural_scope = "cell_type_x_medium",
      fastcore_parallel_task = "cell_type_x_medium",
      parallel_task =
        "directional_reaction_by_matching_metacells_step2",
      vmax_computation_scope = attr(
        vmax_cache, "parallel_scope"
      ) %||% "directional_target_batches_within_shared_models",
      vmax_parallel_tasks = as.integer(
        attr(vmax_cache, "parallel_tasks") %||% length(vmax_cache)
      ),
      vmax_parallel_workers = as.integer(
        attr(vmax_cache, "parallel_workers") %||% 1L
      ),
      vmax_solve_count = length(vmax_cache),
      vmax_reuse_by_cell_type = stats::setNames(
        as.integer(reuse), names(reuse)
      ),
      step2_solver_reuse = paste(
        "one persistent HiGHS model per directional reaction reused across",
        "all matching metacells; one-shot fallback for unsupported backends"
      ),
      worker_cleanup =
        "checkpoint_each_reaction_then_drop_model_and_run_full_gc",
      flux_threshold = flux_threshold,
      scoring_time_limit = "none"
    ),
    method = paste(
      "reaction-parallel microCOMPASS directional LP with persistent native",
      "HiGHS reuse on cell-type-specific medium union GEMs"
    )
  )
}

# Progress-aware entry point; the algorithm remains in the core above.
.rc_run_celltype_microcompass_engine <- function(...) {
  args <- list(...)
  answer <- do.call(
    .rc_run_celltype_microcompass_engine_reaction_core,
    args
  )
  contexts <- .rc_layer2_model_contexts(
    answer$shared_model_cache, mode = "meta_module_gem"
  )
  parts_dir <- .rc_layer2_cache_progress_dir(answer$shared_model_cache)
  run_kind <- .rc_layer2_progress_state$run_kind %||% "primary"
  celltype_col <- as.character(args$celltype_col %||% "cell_type")
  unit_celltypes <- if (celltype_col %in% colnames(answer$unit_meta)) {
    as.character(answer$unit_meta[[celltype_col]])
  } else {
    character()
  }
  for (item in contexts) {
    .rc_layer2_task_event(
      item$context, "penalty_step2_complete", 4L, 4L,
      detail = paste0(
        "matching_units=", sum(unit_celltypes == item$context$cell_type)
      ),
      scope = "scoring", run_kind = run_kind,
      status = "complete", parts_dir = parts_dir
    )
  }
  answer
}
