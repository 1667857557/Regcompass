# Directional microCOMPASS scoring on shared structural models.

.rc_load_microcompass_model <- function(entry, mode) {
  if (identical(mode, "meta_module_gem")) {
    return(.rc_read_celltype_union_gem(
      file = entry$file,
      cell_type = entry$cell_type,
      medium_scenario = entry$medium_scenario,
      expected_checksum = entry$file_checksum %||% NA_character_
    ))
  }
  .rc_cache_gem(entry)
}

.rc_microcompass_object_checksum <- function(x) {
  file <- tempfile("regcompass-contract-", fileext = ".rds")
  on.exit(unlink(file, force = TRUE), add = TRUE)
  saveRDS(x, file, version = 2)
  unname(tools::md5sum(file)[[1L]])
}

.rc_microcompass_model_contract <- function(model_cache, mode) {
  if (!is.list(model_cache) || !length(model_cache)) {
    stop("The shared model cache is empty.", call. = FALSE)
  }
  file_key <- vapply(model_cache, function(entry) {
    as.character(entry$file %||% "")
  }, character(1))
  representative <- match(unique(file_key), file_key)
  rows <- lapply(representative, function(i) {
    entry <- model_cache[[i]]
    model <- .rc_load_microcompass_model(entry, mode)
    validated <- rc_validate_gem(model)
    reaction_order <- colnames(validated$S)
    metabolite_order <- rownames(validated$S)
    data.frame(
      model_file = as.character(entry$file %||% ""),
      model_file_checksum = if (
        !is.null(entry$file) && file.exists(entry$file)
      ) {
        unname(tools::md5sum(entry$file)[[1L]])
      } else {
        NA_character_
      },
      medium_scenario = as.character(entry$medium_scenario),
      n_reactions = ncol(validated$S),
      n_metabolites = nrow(validated$S),
      reaction_order_checksum =
        .rc_microcompass_object_checksum(reaction_order),
      metabolite_order_checksum =
        .rc_microcompass_object_checksum(metabolite_order),
      stoichiometry_bounds_checksum =
        .rc_microcompass_object_checksum(list(
          S = validated$S,
          lb = validated$lb,
          ub = validated$ub
        )),
      stringsAsFactors = FALSE
    )
  })
  answer <- do.call(rbind, rows)
  rownames(answer) <- NULL
  answer[order(answer$medium_scenario, answer$model_file), , drop = FALSE]
}

.rc_full_gem_step2_model_payload <- function(
    model_key, row_ids, model_cache, penalties, vmax_cache,
    omega, solver, flux_threshold, payload_dir) {
  row_ids <- as.character(row_ids)
  if (!length(row_ids)) {
    stop("A full-GEM Step 2 payload requires at least one target row.",
         call. = FALSE)
  }
  first_entry <- model_cache[[row_ids[[1L]]]]
  if (is.null(first_entry)) {
    stop("A full-GEM Step 2 payload references an unknown target row.",
         call. = FALSE)
  }
  model <- .rc_load_microcompass_model(first_entry, "full_gem")
  reactions <- colnames(model$S)
  if (is.null(reactions) || anyNA(reactions) || any(!nzchar(reactions)) ||
      anyDuplicated(reactions)) {
    stop("A full-GEM Step 2 payload has invalid reaction identifiers.",
         call. = FALSE)
  }
  units <- colnames(penalties$penalty)
  if (is.null(units) || anyNA(units) || any(!nzchar(units)) ||
      anyDuplicated(units)) {
    stop("A full-GEM Step 2 payload has invalid unit identifiers.",
         call. = FALSE)
  }
  missing_penalty_rows <- setdiff(reactions, rownames(penalties$penalty))
  if (length(missing_penalty_rows)) {
    stop(
      "Full-GEM Step 2 penalties are missing model reactions: ",
      paste(utils::head(missing_penalty_rows, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  entries <- lapply(row_ids, function(row_id) {
    entry <- model_cache[[row_id]]
    if (is.null(entry) ||
        !identical(as.character(entry$file), as.character(model_key))) {
      stop("A full-GEM Step 2 payload mixes different model files.",
           call. = FALSE)
    }
    list(
      reaction_id = as.character(entry$reaction_id),
      target_direction = as.character(entry$target_direction),
      medium_scenario = as.character(entry$medium_scenario),
      condition = as.character(entry$condition %||% "all")
    )
  })
  names(entries) <- row_ids
  vmax_values <- lapply(row_ids, function(row_id) {
    value <- vmax_cache[[row_id]]
    list(
      feasible = isTRUE(value$feasible),
      vmax = as.numeric(value$vmax),
      status = as.character(value$status),
      flux = numeric()
    )
  })
  names(vmax_values) <- row_ids
  model_checksum <- if (!is.null(first_entry$file) &&
                        file.exists(first_entry$file)) {
    unname(tools::md5sum(first_entry$file)[[1L]])
  } else {
    NA_character_
  }
  payload <- list(
    schema_version = "regcompass_full_gem_step2_compact_payload_v1",
    model = list(
      S = .rc_as_dgCMatrix(model$S),
      lb = model$lb,
      ub = model$ub,
      target_status = model$target_status %||% "not_prechecked",
      file_checksum = model_checksum,
      medium_scenario = as.character(first_entry$medium_scenario),
      condition = as.character(first_entry$condition %||% "all")
    ),
    reactions = reactions,
    units = units,
    penalty = penalties$penalty[reactions, units, drop = FALSE],
    entries = entries,
    vmax = vmax_values,
    omega = as.numeric(omega),
    solver = as.character(solver),
    flux_threshold = as.numeric(flux_threshold)
  )
  token <- substr(.rc_microcompass_object_checksum(list(
    file = as.character(first_entry$file %||% model_key),
    checksum = model_checksum,
    units = units,
    row_ids = row_ids
  )), 1L, 24L)
  file <- file.path(payload_dir, paste0("payload__", token, ".rds"))
  .rc_atomic_save_rds(payload, file)
  rm(model, payload, entries, vmax_values)
  invisible(gc(verbose = FALSE, full = FALSE))
  file
}

.rc_full_gem_step2_reaction_batch_worker <- function(task) {
  if (!is.list(task) ||
      !all(c("payload_file", "row_ids", "checkpoint_dir") %in% names(task))) {
    stop("Malformed full-GEM Step 2 reaction-batch task.", call. = FALSE)
  }
  payload <- readRDS(task$payload_file)
  if (!is.list(payload) ||
      !identical(
        payload$schema_version,
        "regcompass_full_gem_step2_compact_payload_v1"
      )) {
    stop("Malformed full-GEM Step 2 compact payload.", call. = FALSE)
  }
  row_ids <- as.character(task$row_ids)
  if (!length(row_ids) ||
      !all(row_ids %in% names(payload$entries)) ||
      !all(row_ids %in% names(payload$vmax))) {
    stop("Full-GEM Step 2 reaction-batch rows are absent from the payload.",
         call. = FALSE)
  }
  model <- payload$model
  if (!is.list(model) || is.null(model$S) || is.null(model$lb) ||
      is.null(model$ub)) {
    stop("Full-GEM Step 2 compact payload lacks the required LP model state.",
         call. = FALSE)
  }
  if (!identical(colnames(model$S), as.character(payload$reactions))) {
    stop("Full-GEM Step 2 payload reaction order differs from its model.",
         call. = FALSE)
  }
  if (!identical(colnames(payload$penalty), as.character(payload$units)) ||
      !identical(rownames(payload$penalty), as.character(payload$reactions))) {
    stop("Full-GEM Step 2 compact payload penalties are not aligned.",
         call. = FALSE)
  }

  checkpoint_files <- character(length(row_ids))
  step2_engine <- NULL
  on.exit({
    .rc_compass_step2_release_engine(step2_engine)
    rm(model, payload)
    invisible(gc(verbose = FALSE, full = TRUE))
  }, add = TRUE)

  units <- as.character(payload$units)
  for (j in seq_along(row_ids)) {
    row_id <- row_ids[[j]]
    entry <- payload$entries[[row_id]]
    target_index <- match(entry$reaction_id, payload$reactions)
    if (is.na(target_index)) {
      stop("A full-GEM target reaction is absent from its shared model.",
           call. = FALSE)
    }
    prepared <- .rc_compass_step2_prepare(
      S = model$S,
      lb = model$lb,
      ub = model$ub,
      target_reaction = entry$reaction_id,
      vmax_result = payload$vmax[[row_id]],
      target_direction = entry$target_direction,
      omega = payload$omega,
      flux_threshold = payload$flux_threshold
    )
    step2_engine <- if (isTRUE(prepared$runnable)) {
      .rc_compass_step2_new_engine(prepared$template, payload$solver)
    } else {
      NULL
    }

    task_penalty <- rep(NA_real_, length(units))
    task_vmax <- rep(NA_real_, length(units))
    task_feasible <- task_evaluated <- rep(FALSE, length(units))
    names(task_penalty) <- names(task_vmax) <-
      names(task_feasible) <- names(task_evaluated) <- units
    diagnostics <- vector("list", length(units))

    for (i in seq_along(units)) {
      one_unit <- units[[i]]
      unit_penalty <- payload$penalty[, i]
      target_penalty <- unit_penalty[[target_index]]
      evidence_available <- is.finite(unit_penalty)
      solver_penalty <- unit_penalty
      solver_penalty[!evidence_available] <- 0

      if (isTRUE(prepared$runnable)) {
        solved <- .rc_compass_step2_engine_solve(
          step2_engine, solver_penalty
        )
        step2_engine <- solved$engine
        fit <- .rc_compass_step2_result(prepared$template, solved$answer)
      } else {
        fit <- prepared$result
        solved <- NULL
      }

      target_available <- is.finite(target_penalty)
      task_penalty[[one_unit]] <- if (target_available) {
        fit$penalty
      } else {
        NA_real_
      }
      task_vmax[[one_unit]] <- fit$vmax
      task_feasible[[one_unit]] <- isTRUE(fit$feasible)
      task_evaluated[[one_unit]] <- isTRUE(fit$feasible) && target_available
      diagnostics[[i]] <- data.frame(
        row_id = row_id,
        unit_id = one_unit,
        module_id = NA_character_,
        reaction_id = entry$reaction_id,
        target_direction = entry$target_direction,
        medium_scenario = entry$medium_scenario,
        condition = entry$condition,
        strict_feasible = isTRUE(fit$feasible),
        solver_status = fit$solver_status,
        solver_backend = fit$solver_backend %||% "unknown",
        step1_status = fit$step1_status,
        step2_status = fit$step2_status,
        target_status = if (isTRUE(fit$feasible)) {
          "ok"
        } else {
          "medium_directionally_infeasible"
        },
        objective_value = if (target_available) fit$penalty else NA_real_,
        vmax = fit$vmax,
        vmax_reused_from_shared_cache = TRUE,
        step2_model_reused_across_metacells = TRUE,
        target_expression_available = target_available,
        objective_evidence_fraction = mean(evidence_available),
        unavailable_objective_terms = sum(!evidence_available),
        parallel_task = "directional_reaction_x_all_metacells",
        stringsAsFactors = FALSE
      )
      rm(unit_penalty, solver_penalty, fit, solved, evidence_available)
    }

    engine_metrics <- .rc_compass_step2_engine_metrics(step2_engine)
    token <- substr(.rc_microcompass_object_checksum(list(
      row_id = row_id,
      file_checksum = model$file_checksum,
      units = units,
      omega = payload$omega,
      solver = payload$solver,
      flux_threshold = payload$flux_threshold
    )), 1L, 24L)
    checkpoint <- file.path(
      task$checkpoint_dir, paste0("step2__", token, ".rds")
    )
    .rc_atomic_save_rds(list(
      row_id = row_id,
      units = units,
      penalty = task_penalty,
      vmax = task_vmax,
      feasible = task_feasible,
      evaluated = task_evaluated,
      diagnostics = .rc_bind_frames_fill(diagnostics),
      engine_metrics = engine_metrics
    ), checkpoint)
    checkpoint_files[[j]] <- checkpoint

    .rc_compass_step2_release_engine(step2_engine)
    step2_engine <- NULL
    rm(
      diagnostics, task_penalty, task_vmax, task_feasible,
      task_evaluated, prepared, engine_metrics
    )
    invisible(gc(verbose = FALSE, full = FALSE))
  }
  checkpoint_files
}

.rc_run_shared_full_gem_engine_core <- function(
    layer1, gem, target_reactions = NULL,
    medium_table = NULL, medium_scenarios = NULL,
    mode = c("full_gem", "meta_module_gem"),
    reaction_membership = NULL, core_reactions = NULL,
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
  mode <- match.arg(mode)
  if (!identical(mode, "full_gem")) {
    stop("The shared full-GEM engine accepts only `full_gem` mode.",
         call. = FALSE)
  }
  .rc_validate_full_gem_model_params(model_params)
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
  gem <- rc_annotate_reaction_roles(gem)
  direction_diagnostics <- NULL

  if (!is.null(model_cache_override)) {
    if (!is.list(model_cache_override) || !length(model_cache_override) ||
        is.null(attr(model_cache_override, "summary")) ||
        !identical(attr(model_cache_override, "completion_method"), "none") ||
        !identical(attr(model_cache_override, "fastcore_executed"), FALSE) ||
        !identical(attr(model_cache_override, "corda2_executed"), FALSE) ||
        !identical(
          attr(model_cache_override, "medium_handling"),
          "exchange_bounds_only_no_reaction_deletion"
        )) {
      stop(
        "`model_cache_override` is not an audited COMPASS-style ",
        "medium-constrained full-GEM cache.",
        call. = FALSE
      )
    }
    model_cache <- model_cache_override
    directions <- unique(do.call(rbind, lapply(
      model_cache,
      function(entry) {
        data.frame(
          reaction_id = as.character(entry$reaction_id),
          target_direction = as.character(entry$target_direction),
          medium_scenario = as.character(entry$medium_scenario),
          stringsAsFactors = FALSE
        )
      }
    )))
    cached_media <- unique(as.character(directions$medium_scenario))
    requested_media <- unique(as.character(
      medium_scenarios$medium_scenario_id
    ))
    if (!setequal(cached_media, requested_media)) {
      stop("The shared cache and requested media differ.", call. = FALSE)
    }
  } else {
    if (is.null(target_reactions) || !length(target_reactions)) {
      stop("`target_reactions` is required in full-GEM mode.", call. = FALSE)
    }
    directions <- rc_prepare_directional_targets(
      gem, target_reactions, target_direction
    )
    direction_diagnostics <- directions
    directions <- directions[
      directions$target_direction %in% c("forward", "reverse"),
      , drop = FALSE
    ]
    if (!nrow(directions)) {
      stop("No target reaction directions are allowed by the GEM bounds.",
           call. = FALSE)
    }
    model_cache <- rc_build_full_gem_cache(
      gem = gem,
      dirs = directions,
      medium_scenarios = medium_scenarios,
      cache_dir = model_params$cache_dir %||%
        tempfile("RegCompassR_full_gem_cache_"),
      conditions = "all",
      solver = solver,
      time_limit = model_params$completion_time_limit %||% NA_real_,
      flux_consistency_epsilon = flux_threshold
    )
    directions <- unique(do.call(rbind, lapply(
      model_cache,
      function(entry) {
        data.frame(
          reaction_id = as.character(entry$reaction_id),
          target_direction = as.character(entry$target_direction),
          medium_scenario = as.character(entry$medium_scenario),
          stringsAsFactors = FALSE
        )
      }
    )))
  }

  units <- colnames(matrices$reaction_expression)
  unit_meta <- matrices$unit_meta
  row_ids <- names(model_cache)
  model_keys <- vapply(row_ids, function(row_id) {
    entry <- model_cache[[row_id]]
    entry$file %||% paste0("memory::", row_id)
  }, character(1))
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
    colnames(.rc_load_microcompass_model(entry, mode)$S)
  }), use.names = FALSE))
  penalties <- rc_compute_multiome_penalty(
    rc_align_reaction_expression(
      matrices$reaction_expression, all_reactions, NA_real_
    ),
    reaction_roles = gem$reaction_roles
  )
  vmax_cache <- .rc_build_microcompass_vmax_cache(
    model_cache = model_cache,
    mode = mode,
    model_keys = model_keys,
    solver = solver,
    flux_threshold = flux_threshold,
    parallel = parallel,
    BPPARAM = BPPARAM
  )
  vmax_parallel_tasks <- as.integer(
    attr(vmax_cache, "parallel_tasks") %||% length(vmax_cache)
  )
  vmax_parallel_workers <- as.integer(
    attr(vmax_cache, "parallel_workers") %||% 1L
  )
  vmax_computation_scope <- attr(vmax_cache, "parallel_scope") %||%
    "shared_model_x_directional_target_once"
  vmax_solve_count <- length(vmax_cache)
  vmax_cache_diagnostics <- do.call(rbind, lapply(row_ids, function(row_id) {
    entry <- model_cache[[row_id]]
    value <- vmax_cache[[row_id]]
    data.frame(
      row_id = row_id,
      reaction_id = as.character(entry$reaction_id),
      target_direction = as.character(entry$target_direction),
      medium_scenario = as.character(entry$medium_scenario),
      vmax = as.numeric(value$vmax),
      feasible = isTRUE(value$feasible),
      status = as.character(value$status),
      computation_scope = "shared_model_x_directional_target_once",
      stringsAsFactors = FALSE
    )
  }))
  rownames(vmax_cache_diagnostics) <- NULL

  penalty <- vmax <- matrix(
    NA_real_,
    nrow = length(row_ids),
    ncol = length(units),
    dimnames = list(row_ids, units)
  )
  feasible <- evaluated <- matrix(
    FALSE,
    nrow = length(row_ids),
    ncol = length(units),
    dimnames = list(row_ids, units)
  )

  checkpoint_root <- model_params$cache_dir %||%
    dirname(model_cache[[row_ids[[1L]]]]$file)
  dir.create(checkpoint_root, recursive = TRUE, showWarnings = FALSE)
  checkpoint_dir <- tempfile(
    pattern = "full_gem_step2_reaction_tasks_",
    tmpdir = checkpoint_root
  )
  payload_dir <- file.path(checkpoint_dir, "compact_payloads")
  dir.create(payload_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(checkpoint_dir, recursive = TRUE, force = TRUE), add = TRUE)

  step2_workers <- .rc_microcompass_worker_count(
    parallel = parallel,
    BPPARAM = BPPARAM,
    n_tasks = length(row_ids)
  )
  batch_specs <- .rc_step2_model_batches(model_keys, step2_workers)
  step2_task_count <- length(batch_specs)
  payload_files <- vapply(unique_model_keys, function(model_key) {
    model_rows <- names(model_keys)[model_keys == model_key]
    .rc_full_gem_step2_model_payload(
      model_key = model_key,
      row_ids = model_rows,
      model_cache = model_cache,
      penalties = penalties,
      vmax_cache = vmax_cache,
      omega = omega,
      solver = solver,
      flux_threshold = flux_threshold,
      payload_dir = payload_dir
    )
  }, character(1))
  names(payload_files) <- unique_model_keys
  tasks <- lapply(batch_specs, function(spec) {
    list(
      payload_file = payload_files[[as.character(spec$model_key)]],
      row_ids = as.character(spec$row_ids),
      checkpoint_dir = checkpoint_dir
    )
  })
  names(tasks) <- names(batch_specs)

  penalties$penalty <- NULL
  rm(vmax_cache, matrices, all_reactions)
  invisible(gc(verbose = FALSE, full = TRUE))

  checkpoint_groups <- rc_parallel_lapply(
    tasks,
    .rc_full_gem_step2_reaction_batch_worker,
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )
  checkpoint_files <- unlist(checkpoint_groups, use.names = FALSE)
  if (length(checkpoint_files) != length(row_ids)) {
    stop("Full-GEM Step 2 reaction batches returned an incomplete checkpoint set.",
         call. = FALSE)
  }

  diagnostics <- vector("list", length(checkpoint_files))
  step2_engine_metrics <- vector("list", length(checkpoint_files))
  observed_rows <- character(length(checkpoint_files))
  for (i in seq_along(checkpoint_files)) {
    result <- readRDS(checkpoint_files[[i]])
    row_id <- as.character(result$row_id)
    if (!row_id %in% row_ids ||
        !identical(names(result$penalty), as.character(result$units))) {
      stop("A full-GEM reaction-level Step 2 checkpoint is malformed.",
           call. = FALSE)
    }
    observed_rows[[i]] <- row_id
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
    rm(result, metrics)
    unlink(checkpoint_files[[i]], force = TRUE)
    invisible(gc(verbose = FALSE, full = FALSE))
  }
  if (anyDuplicated(observed_rows) || !setequal(observed_rows, row_ids)) {
    stop("Full-GEM Step 2 reaction batches did not score every target exactly once.",
         call. = FALSE)
  }
  unlink(payload_files, force = TRUE)
  rm(checkpoint_groups, checkpoint_files, tasks, batch_specs)
  invisible(gc(verbose = FALSE, full = TRUE))

  score <- rc_compass_score_from_penalty(penalty, feasible)
  lp_diagnostics <- .rc_bind_frames_fill(diagnostics)
  model_diagnostics <- .rc_bind_frames_fill(lapply(
    representative_rows,
    function(row_id) {
      entry <- model_cache[[row_id]]
      model <- .rc_load_microcompass_model(entry, mode)
      out <- model$closure_diagnostics %||% data.frame()
      rm(model)
      out
    }
  ))

  list(
    score = score,
    penalty = penalty,
    vmax = vmax,
    feasible = feasible,
    evaluated = evaluated,
    target_direction = directions,
    direction_diagnostics = direction_diagnostics,
    medium_scenarios = medium_scenarios,
    model_mode = mode,
    shared_model_cache = model_cache,
    model_cache_summary = attr(model_cache, "summary"),
    structural_model_contract =
      .rc_microcompass_model_contract(model_cache, mode),
    model_diagnostics = model_diagnostics,
    vmax_cache_diagnostics = vmax_cache_diagnostics,
    lp_diagnostics = lp_diagnostics,
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
      shared_gem_scope = paste(
        "one complete medium-constrained full GEM per medium",
        "shared across all units"
      ),
      structural_scope = "medium_x_complete_full_gem",
      model_completion = "none",
      structural_completion = "none",
      structural_completion_algorithm = "compass_medium_bounds_only",
      medium_handling = "exchange_bounds_only_no_reaction_deletion",
      medium_direct_reaction_deletion = FALSE,
      target_feasibility = paste(
        "directional vmax under the medium; skip Step 2 when",
        "vmax is below flux_threshold"
      ),
      fastcore_executed = FALSE,
      corda2_executed = FALSE,
      parallel_task = "directional_reaction_by_all_metacells_step2",
      step2_dispatch = "model_scoped_reaction_batches",
      step2_parallel_tasks = as.integer(step2_task_count),
      step2_parallel_workers = as.integer(step2_workers),
      step2_worker_payload = paste(
        "file-backed compact S/lb/ub, full-GEM penalties, target metadata",
        "and cached vmax; workers return checkpoint paths only"
      ),
      step2_model_load_reuse = paste(
        "each reaction batch loads one compact full-GEM payload once and",
        "reuses it across its directional reactions"
      ),
      step2_solver_reuse = paste(
        "one persistent HiGHS model per directional reaction reused across",
        "all metacells by objective-coefficient updates"
      ),
      worker_cleanup = paste(
        "checkpoint each reaction, release its HiGHS engine, remove",
        "reaction-local objects, then release the dispatch-local worker pool"
      ),
      vmax_computation_scope = vmax_computation_scope,
      vmax_parallel_tasks = vmax_parallel_tasks,
      vmax_parallel_workers = vmax_parallel_workers,
      vmax_solve_count = vmax_solve_count,
      vmax_reuse_factor = length(units),
      flux_threshold = flux_threshold,
      scoring_time_limit = "none"
    ),
    method = paste(
      "microCOMPASS shared medium-constrained full-GEM directional LP with",
      "reaction-parallel persistent HiGHS Step 2"
    )
  )
}

.rc_run_microcompass_engine <- function(
    layer1, gem, target_reactions = NULL,
    medium_table = NULL, medium_scenarios = NULL,
    mode = c("full_gem", "meta_module_gem"),
    reaction_membership = NULL, core_reactions = NULL,
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
  mode <- match.arg(mode)
  unit <- match.arg(unit)
  solver <- match.arg(solver)
  target_direction <- match.arg(target_direction)
  if (identical(mode, "meta_module_gem")) {
    return(.rc_run_celltype_microcompass_engine(
      layer1 = layer1,
      gem = gem,
      target_reactions = target_reactions,
      medium_table = medium_table,
      medium_scenarios = medium_scenarios,
      reaction_membership = reaction_membership,
      core_reactions = core_reactions,
      unit = unit,
      condition_col = condition_col,
      sample_col = sample_col,
      celltype_col = celltype_col,
      model_params = model_params,
      omega = omega,
      target_direction = target_direction,
      parallel = parallel,
      solver = solver,
      flux_threshold = flux_threshold,
      BPPARAM = BPPARAM,
      model_cache_override = model_cache_override
    ))
  }
  .rc_run_shared_full_gem_engine(
    layer1 = layer1,
    gem = gem,
    target_reactions = target_reactions,
    medium_table = medium_table,
    medium_scenarios = medium_scenarios,
    mode = "full_gem",
    reaction_membership = NULL,
    core_reactions = NULL,
    unit = unit,
    condition_col = condition_col,
    sample_col = sample_col,
    celltype_col = celltype_col,
    model_params = model_params,
    omega = omega,
    target_direction = target_direction,
    parallel = parallel,
    solver = solver,
    flux_threshold = flux_threshold,
    BPPARAM = BPPARAM,
    model_cache_override = model_cache_override
  )
}

#' Run directional minimum-evidence-discordance LPs
rc_run_microcompass <- function(
    layer1, gem, target_reactions = NULL,
    medium_table = NULL, medium_scenarios = NULL,
    mode = c("full_gem", "meta_module_gem"),
    reaction_membership = NULL, core_reactions = NULL,
    unit = c("metacell", "sample_celltype"),
    condition_col = "condition", sample_col = NULL,
    celltype_col = "cell_type", model_params = list(),
    omega = 0.95,
    target_direction = c("both", "forward", "reverse"),
    parallel = TRUE,
    solver = c("highs", "gurobi", "glpk"),
    flux_threshold = 1e-8,
    BPPARAM = NULL) {
  mode <- match.arg(mode)
  unit <- match.arg(unit)
  target_direction <- match.arg(target_direction)
  solver <- match.arg(solver)
  if (!is.list(model_params)) {
    stop("`model_params` must be a list.", call. = FALSE)
  }
  if (identical(unit, "metacell")) {
    warning(
      paste(
        "Metacells are valid within-dataset statistical units, but their",
        "P values are not sample-level biological-replicate inference."
      ),
      call. = FALSE
    )
  }
  if (is.null(medium_scenarios) && is.null(medium_table)) {
    model_info <- gem$model_info %||% list()
    recorded_species <- tolower(as.character(model_info$species %||% ""))
    recorded_source <- tolower(as.character(model_info$source %||% ""))
    species_provenance <- recorded_species %in% c("human", "mouse") ||
      grepl("human-gem|mouse-gem", recorded_source)
    if (isTRUE(species_provenance)) {
      species <- if (
        identical(recorded_species, "mouse") ||
          grepl("mouse-gem", recorded_source)
      ) {
        "mouse"
      } else {
        "human"
      }
      medium_scenarios <- rc_make_medium_scenarios(
        gem,
        scenario = if (identical(species, "mouse")) {
          "mouse_plasma"
        } else {
          "normal_human_plasma"
        },
        species = species,
        strict_preset_matching = TRUE
      )
    } else {
      stop(
        "GEM species provenance is unavailable; supply an explicit, ",
        "biologically justified `medium_scenarios` or `medium_table`.",
        call. = FALSE
      )
    }
  }
  cache_dir <- model_params$cache_dir %||% NULL
  if (!is.null(cache_dir)) {
    if (!is.character(cache_dir) || length(cache_dir) != 1L ||
        is.na(cache_dir) || !nzchar(cache_dir)) {
      stop("`model_params$cache_dir` must be one non-empty path.",
           call. = FALSE)
    }
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  answer <- .rc_run_microcompass_engine(
    layer1 = layer1,
    gem = gem,
    target_reactions = target_reactions,
    medium_table = medium_table,
    medium_scenarios = medium_scenarios,
    mode = mode,
    reaction_membership = reaction_membership,
    core_reactions = core_reactions,
    unit = unit,
    condition_col = condition_col,
    sample_col = sample_col,
    celltype_col = celltype_col,
    model_params = model_params,
    omega = omega,
    target_direction = target_direction,
    parallel = parallel,
    solver = solver,
    flux_threshold = flux_threshold,
    BPPARAM = BPPARAM
  )
  answer$relative_penalty_rank <- answer$score
  answer$score_semantics <- attr(answer$score, "score_semantics") %||%
    "within_target_relative_penalty_rank_not_probability"
  answer$noninformative_target <- attr(
    answer$score, "noninformative_target"
  )
  answer$primary_output <- "penalty"
  answer$primary_output_semantics <-
    "minimum evidence-discordance penalty; lower means stronger support"
  answer$params$inference_unit <- unit
  answer$params$model_cache_dir <- cache_dir
  if (!is.null(cache_dir)) {
    saveRDS(
      answer$model_cache_summary,
      file.path(cache_dir, "model_cache_summary.rds")
    )
    saveRDS(
      answer$model_diagnostics,
      file.path(cache_dir, "model_diagnostics.rds")
    )
    saveRDS(
      answer$vmax_cache_diagnostics,
      file.path(cache_dir, "vmax_cache_diagnostics.rds")
    )
    saveRDS(
      answer$lp_diagnostics,
      file.path(cache_dir, "lp_diagnostics.rds")
    )
    saveRDS(
      answer$medium_scenarios,
      file.path(cache_dir, "medium_scenarios.rds")
    )
    saveRDS(
      answer$params,
      file.path(cache_dir, "microcompass_parameters.rds")
    )
    files <- if (is.data.frame(answer$model_cache_summary) &&
                 "file" %in% colnames(answer$model_cache_summary)) {
      unique(as.character(answer$model_cache_summary$file))
    } else {
      character()
    }
    files <- files[file.exists(files)]
    manifest <- data.frame(
      file = files,
      size_bytes = as.numeric(file.info(files)$size),
      checksum = if (length(files)) {
        unname(tools::md5sum(files))
      } else {
        character()
      },
      stringsAsFactors = FALSE
    )
    saveRDS(manifest, file.path(cache_dir, "model_file_manifest.rds"))
    answer$model_file_manifest <- manifest
  }
  answer
}

#' Summarize a microCOMPASS result

# Progress-aware entry point; the algorithm remains in the core above.
.rc_run_shared_full_gem_engine <- function(...) {
  progress_state <- get0(
    ".rc_layer2_progress_state", mode = "environment", inherits = TRUE
  )
  answer <- do.call(
    .rc_run_shared_full_gem_engine_core,
    list(...)
  )
  contexts <- .rc_layer2_model_contexts(
    answer$shared_model_cache, mode = "full_gem"
  )
  parts_dir <- .rc_layer2_cache_progress_dir(answer$shared_model_cache)
  run_kind <- progress_state$run_kind %||% "primary"
  for (item in contexts) {
    .rc_layer2_task_event(
      item$context, "penalty_step2_complete", 4L, 4L,
      detail = paste0("units=", nrow(answer$unit_meta)),
      scope = "scoring", run_kind = run_kind,
      status = "complete", parts_dir = parts_dir
    )
  }
  answer
}
