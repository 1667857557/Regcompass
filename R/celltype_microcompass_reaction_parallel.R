# Reaction-granular directional LP scoring on cell-type-specific models.

.rc_step2_compact_vmax_value <- function(value) {
  list(
    feasible = isTRUE(value$feasible),
    vmax = as.numeric(value$vmax),
    status = as.character(value$status),
    flux = numeric()
  )
}

.rc_step2_model_batches <- function(model_keys, workers) {
  if (is.null(names(model_keys)) || any(!nzchar(names(model_keys)))) {
    stop("Step 2 model keys require named directional target rows.",
         call. = FALSE)
  }
  unique_keys <- unique(unname(model_keys))
  n_models <- length(unique_keys)
  workers <- max(1L, as.integer(workers[[1L]]))
  rows_by_model <- lapply(unique_keys, function(model_key) {
    names(model_keys)[model_keys == model_key]
  })
  n_rows <- vapply(rows_by_model, length, integer(1))
  target_batches <- min(length(model_keys), max(n_models, workers))
  n_batches <- rep.int(1L, n_models)
  remaining <- target_batches - n_models
  while (remaining > 0L) {
    eligible <- which(n_batches < n_rows)
    if (!length(eligible)) break
    load <- n_rows[eligible] / n_batches[eligible]
    chosen <- eligible[[which.max(load)]]
    n_batches[[chosen]] <- n_batches[[chosen]] + 1L
    remaining <- remaining - 1L
  }

  tasks <- list()
  cursor <- 0L
  for (model_index in seq_along(unique_keys)) {
    selected_rows <- rows_by_model[[model_index]]
    batches <- n_batches[[model_index]]
    batch_id <- ceiling(
      seq_along(selected_rows) * batches / length(selected_rows)
    )
    groups <- split(selected_rows, batch_id)
    for (batch_index in seq_along(groups)) {
      cursor <- cursor + 1L
      tasks[[cursor]] <- list(
        model_key = unique_keys[[model_index]],
        row_ids = as.character(groups[[batch_index]])
      )
      names(tasks)[[cursor]] <- paste0(
        "model_", model_index, "__batch_", batch_index
      )
    }
  }
  tasks
}

.rc_step2_model_payload <- function(
    model_key, row_ids, model_cache, unit_celltype, penalties, vmax_cache,
    omega, solver, flux_threshold, payload_dir,
    control_penalties = NULL) {
  row_ids <- as.character(row_ids)
  first_entry <- model_cache[[row_ids[[1L]]]]
  eligible <- names(unit_celltype)[unit_celltype == first_entry$cell_type]
  if (!length(eligible)) {
    stop("No Layer 1 units match cell type `", first_entry$cell_type, "`.",
         call. = FALSE)
  }
  model <- .rc_read_celltype_union_gem(
    first_entry$file, first_entry$cell_type,
    first_entry$medium_scenario, first_entry$file_checksum
  )
  reactions <- colnames(model$S)
  if (is.null(reactions) || anyNA(reactions) || any(!nzchar(reactions))) {
    stop("A cell-type union GEM has invalid reaction identifiers.",
         call. = FALSE)
  }
  entries <- lapply(row_ids, function(row_id) {
    entry <- model_cache[[row_id]]
    if (!identical(as.character(entry$file), as.character(model_key))) {
      stop("A Step 2 payload mixes different union GEM files.", call. = FALSE)
    }
    list(
      reaction_id = as.character(entry$reaction_id),
      target_direction = as.character(entry$target_direction),
      cell_type = as.character(entry$cell_type),
      medium_scenario = as.character(entry$medium_scenario)
    )
  })
  names(entries) <- row_ids
  vmax_values <- lapply(row_ids, function(row_id) {
    .rc_step2_compact_vmax_value(vmax_cache[[row_id]])
  })
  names(vmax_values) <- row_ids
  penalty_matrix <- penalties$penalty[reactions, eligible, drop = FALSE]
  penalty_evidence <- .rc_step2_penalty_evidence_stats(penalty_matrix)
  control_penalty_matrix <- NULL
  control_penalty_evidence <- NULL
  control_identical <- FALSE
  if (!is.null(control_penalties)) {
    if (!identical(colnames(control_penalties$penalty), colnames(penalties$penalty)) ||
        !all(reactions %in% rownames(control_penalties$penalty))) {
      stop("RNA-control penalties are not aligned to primary penalties.",
           call. = FALSE)
    }
    control_penalty_matrix <-
      control_penalties$penalty[reactions, eligible, drop = FALSE]
    control_penalty_evidence <-
      .rc_step2_penalty_evidence_stats(control_penalty_matrix)
    control_identical <- vapply(
      seq_len(ncol(penalty_matrix)),
      function(i) identical(
        control_penalty_matrix[, i, drop = TRUE],
        penalty_matrix[, i, drop = TRUE]
      ),
      logical(1)
    )
    names(control_identical) <- colnames(penalty_matrix)
  }
  payload <- list(
    schema_version = "regcompass_step2_compact_payload_v1",
    model = list(
      S = model$S,
      lb = model$lb,
      ub = model$ub,
      target_status = model$target_status %||% NA_character_,
      file_checksum = as.character(first_entry$file_checksum),
      cell_type = as.character(first_entry$cell_type),
      medium_scenario = as.character(first_entry$medium_scenario)
    ),
    reactions = reactions,
    units = eligible,
    penalty = penalty_matrix,
    penalty_evidence = penalty_evidence,
    control_penalty = control_penalty_matrix,
    control_penalty_evidence = control_penalty_evidence,
    control_identical = control_identical,
    entries = entries,
    vmax = vmax_values,
    omega = as.numeric(omega),
    solver = as.character(solver),
    flux_threshold = as.numeric(flux_threshold)
  )
  token <- substr(.rc_microcompass_object_checksum(list(
    file = first_entry$file,
    checksum = first_entry$file_checksum,
    units = eligible,
    row_ids = row_ids
  )), 1L, 24L)
  file <- file.path(payload_dir, paste0("payload__", token, ".rds"))
  .rc_atomic_save_rds(payload, file)
  rm(
    model, payload, entries, vmax_values, penalty_matrix, penalty_evidence,
    control_penalty_matrix, control_penalty_evidence
  )
  invisible(gc(verbose = FALSE, full = FALSE))
  file
}

.rc_step2_reaction_batch_worker <- function(task) {
  if (!is.list(task) ||
      !all(c("payload_file", "row_ids", "checkpoint_dir") %in% names(task))) {
    stop("Malformed Step 2 reaction-batch task.", call. = FALSE)
  }
  payload <- readRDS(task$payload_file)
  if (!is.list(payload) ||
      !identical(
        payload$schema_version,
        "regcompass_step2_compact_payload_v1"
      )) {
    stop("Malformed Step 2 compact payload.", call. = FALSE)
  }
  row_ids <- as.character(task$row_ids)
  if (!length(row_ids) ||
      !all(row_ids %in% names(payload$entries)) ||
      !all(row_ids %in% names(payload$vmax))) {
    stop("Step 2 reaction-batch rows are absent from the compact payload.",
         call. = FALSE)
  }
  model <- payload$model
  if (!is.list(model) || is.null(model$S) || is.null(model$lb) ||
      is.null(model$ub)) {
    stop("Step 2 compact payload lacks the required LP model state.",
         call. = FALSE)
  }
  if (!identical(colnames(model$S), as.character(payload$reactions))) {
    stop("Step 2 compact payload reaction order differs from its union GEM.",
         call. = FALSE)
  }
  if (!identical(colnames(payload$penalty), as.character(payload$units)) ||
      !identical(rownames(payload$penalty), as.character(payload$reactions))) {
    stop("Step 2 compact payload penalties are not aligned.",
         call. = FALSE)
  }
  paired_control <- !is.null(payload$control_penalty)
  if (paired_control &&
      (!identical(dimnames(payload$control_penalty), dimnames(payload$penalty)))) {
    stop("Paired RNA-control penalties are not aligned.", call. = FALSE)
  }

  checkpoint_files <- character(length(row_ids))
  step2_engine <- NULL
  on.exit({
    .rc_compass_step2_release_engine(step2_engine)
    rm(model, payload)
    invisible(gc(verbose = FALSE, full = TRUE))
  }, add = TRUE)

  units <- as.character(payload$units)
  evidence <- payload$penalty_evidence %||%
    .rc_step2_penalty_evidence_stats(payload$penalty)
  control_evidence <- if (paired_control) {
    payload$control_penalty_evidence %||%
      .rc_step2_penalty_evidence_stats(payload$control_penalty)
  } else {
    NULL
  }
  control_reuse_mask <- if (paired_control) {
    value <- as.logical(payload$control_identical)
    if (length(value) != length(units) || anyNA(value)) {
      stop("Paired RNA-control exact-reuse mask is not aligned to units.",
           call. = FALSE)
    }
    names(value) <- units
    value
  } else {
    setNames(rep(FALSE, length(units)), units)
  }

  for (j in seq_along(row_ids)) {
    row_id <- row_ids[[j]]
    entry <- payload$entries[[row_id]]
    target_index <- match(entry$reaction_id, payload$reactions)
    if (is.na(target_index)) {
      stop("A target reaction is absent from its cell-type union GEM.",
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
    new_engine <- function() {
      if (!isTRUE(prepared$runnable)) return(NULL)
      .rc_compass_step2_new_engine(
        prepared$template, payload$solver,
        persistent_required = identical(payload$solver, "highs")
      )
    }
    step2_engine <- new_engine()
    primary <- .rc_compass_step2_route_solve(
      prepared, step2_engine, payload$penalty, evidence,
      target_index, units
    )
    step2_engine <- primary$engine
    primary_metrics <- .rc_compass_step2_engine_metrics(step2_engine)

    control <- NULL
    control_metrics <- list(
      engine = "not_run", n_solves = 0L,
      n_objective_updates = 0L, n_fallback = 0L
    )
    reused_control <- control_reuse_mask
    if (paired_control) {
      control <- .rc_compass_step2_route_solve(
        prepared, step2_engine, payload$control_penalty,
        control_evidence, target_index, units,
        reuse_mask = reused_control,
        reuse_result = primary
      )
      step2_engine <- control$engine
      control_metrics <- .rc_compass_step2_engine_metrics_delta(
        .rc_compass_step2_engine_metrics(step2_engine),
        primary_metrics
      )
    }

    diagnostics <- data.frame(
      row_id = rep(row_id, length(units)),
      unit_id = units,
      module_id = rep("CELLTYPE_MEDIUM_UNION_GEM", length(units)),
      cell_type = rep(entry$cell_type, length(units)),
      reaction_id = rep(entry$reaction_id, length(units)),
      target_direction = rep(entry$target_direction, length(units)),
      medium_scenario = rep(entry$medium_scenario, length(units)),
      condition = rep("all", length(units)),
      strict_feasible = unname(primary$feasible),
      solver_status = primary$solver_status,
      solver_backend = primary$solver_backend,
      step1_status = primary$step1_status,
      step2_status = primary$step2_status,
      target_status = ifelse(
        primary$feasible, "ok", "structurally_infeasible"
      ),
      objective_value = unname(primary$penalty),
      vmax = unname(primary$vmax),
      vmax_reused_from_celltype_cache = rep(TRUE, length(units)),
      step2_model_reused_across_metacells = rep(TRUE, length(units)),
      target_expression_available = primary$target_available,
      objective_evidence_fraction = evidence$fraction,
      unavailable_objective_terms = evidence$unavailable,
      parallel_task = rep(
        "directional_reaction_x_matching_metacells", length(units)
      ),
      stringsAsFactors = FALSE
    )
    control_diagnostics <- if (paired_control) {
      data.frame(
        row_id = rep(row_id, length(units)),
        unit_id = units,
        cell_type = rep(entry$cell_type, length(units)),
        reaction_id = rep(entry$reaction_id, length(units)),
        target_direction = rep(entry$target_direction, length(units)),
        medium_scenario = rep(entry$medium_scenario, length(units)),
        objective_value = unname(control$penalty),
        strict_feasible = unname(control$feasible),
        solver_status = control$solver_status,
        solver_backend = control$solver_backend,
        objective_identical_to_primary = unname(reused_control),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame()
    }

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
      penalty = primary$penalty,
      vmax = primary$vmax,
      feasible = primary$feasible,
      evaluated = primary$evaluated,
      diagnostics = diagnostics,
      engine_metrics = primary_metrics,
      control = if (paired_control) list(
        penalty = control$penalty,
        vmax = control$vmax,
        feasible = control$feasible,
        evaluated = control$evaluated,
        diagnostics = control_diagnostics,
        engine_metrics = control_metrics,
        reused_from_primary = isTRUE(all(reused_control)),
        reused_from_primary_by_unit = reused_control,
        shared_target_engine = TRUE
      ) else NULL
    ), checkpoint)
    checkpoint_files[[j]] <- checkpoint

    .rc_compass_step2_release_engine(step2_engine)
    step2_engine <- NULL
    rm(
      diagnostics, control_diagnostics, prepared, primary, control,
      primary_metrics, control_metrics
    )
    invisible(gc(verbose = FALSE, full = FALSE))
  }
  checkpoint_files
}

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
    model_cache_override = NULL,
    control_layer1 = NULL) {
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
  control_matrices <- if (!is.null(control_layer1)) {
    value <- rc_layer2_unit_matrices(
      control_layer1,
      if (identical(unit, "metacell")) "metacell" else "sample_celltype",
      sample_col, celltype_col, condition_col
    )
    if (!identical(dimnames(value$reaction_expression),
                   dimnames(matrices$reaction_expression)) ||
        !identical(value$unit_meta, matrices$unit_meta)) {
      stop("Paired RNA-control units do not align with primary Layer 2 units.",
           call. = FALSE)
    }
    value
  } else NULL
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
  control_penalties <- if (!is.null(control_matrices)) {
    rc_compute_multiome_penalty(
      rc_align_reaction_expression(
        control_matrices$reaction_expression, all_reactions, NA_real_
      ),
      reaction_roles = gem$reaction_roles
    )
  } else NULL
  vmax_cache <- .rc_build_microcompass_vmax_cache(
    model_cache = model_cache,
    mode = "meta_module_gem",
    model_keys = model_keys,
    solver = solver,
    flux_threshold = flux_threshold,
    parallel = parallel,
    BPPARAM = BPPARAM
  )
  vmax_computation_scope <- attr(vmax_cache, "parallel_scope") %||%
    "directional_target_batches_within_shared_models"
  vmax_parallel_tasks <- as.integer(
    attr(vmax_cache, "parallel_tasks") %||% length(vmax_cache)
  )
  vmax_parallel_workers <- as.integer(
    attr(vmax_cache, "parallel_workers") %||% 1L
  )
  vmax_solve_count <- length(vmax_cache)
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
  control_penalty <- if (!is.null(control_penalties)) penalty else NULL
  control_feasible <- if (!is.null(control_penalties)) feasible else NULL
  control_evaluated <- if (!is.null(control_penalties)) evaluated else NULL

  checkpoint_root <- model_params$cache_dir %||%
    dirname(model_cache[[row_ids[[1L]]]]$file)
  checkpoint_dir <- tempfile(
    pattern = "step2_reaction_tasks_",
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
    .rc_step2_model_payload(
      model_key = model_key,
      row_ids = model_rows,
      model_cache = model_cache,
      unit_celltype = unit_celltype,
      penalties = penalties,
      vmax_cache = vmax_cache,
      omega = omega,
      solver = solver,
      flux_threshold = flux_threshold,
      payload_dir = payload_dir,
      control_penalties = control_penalties
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
  if (!is.null(control_penalties)) control_penalties$penalty <- NULL
  rm(vmax_cache, matrices, control_matrices, all_reactions)
  invisible(gc(verbose = FALSE, full = TRUE))

  checkpoint_groups <- rc_parallel_lapply(
    tasks,
    .rc_step2_reaction_batch_worker,
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )
  checkpoint_files <- unlist(checkpoint_groups, use.names = FALSE)
  if (length(checkpoint_files) != length(row_ids)) {
    stop("Step 2 reaction batches returned an incomplete checkpoint set.",
         call. = FALSE)
  }
  diagnostics <- vector("list", length(checkpoint_files))
  step2_engine_metrics <- vector("list", length(checkpoint_files))
  control_diagnostics <- vector("list", length(checkpoint_files))
  control_step2_engine_metrics <- vector("list", length(checkpoint_files))
  control_reused <- logical(length(checkpoint_files))
  observed_rows <- character(length(checkpoint_files))
  for (i in seq_along(checkpoint_files)) {
    result <- readRDS(checkpoint_files[[i]])
    row_id <- as.character(result$row_id)
    if (!row_id %in% row_ids ||
        !identical(names(result$penalty), as.character(result$units))) {
      stop("A reaction-level Step 2 checkpoint is malformed.",
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
    if (!is.null(control_penalty)) {
      ctrl <- result$control
      if (!is.list(ctrl) ||
          !identical(names(ctrl$penalty), as.character(result$units))) {
        stop("A paired RNA-control checkpoint is malformed.", call. = FALSE)
      }
      control_penalty[row_id, result$units] <- ctrl$penalty
      control_feasible[row_id, result$units] <- ctrl$feasible
      control_evaluated[row_id, result$units] <- ctrl$evaluated
      control_diagnostics[[i]] <- ctrl$diagnostics
      reuse_mask <- as.logical(
        ctrl$reused_from_primary_by_unit %||%
          rep(isTRUE(ctrl$reused_from_primary), length(result$units))
      )
      if (length(reuse_mask) != length(result$units) || anyNA(reuse_mask)) {
        stop("A paired RNA-control checkpoint has a malformed reuse mask.",
             call. = FALSE)
      }
      control_reused[[i]] <- isTRUE(all(reuse_mask))
      control_metrics <- ctrl$engine_metrics %||% list()
      control_step2_engine_metrics[[i]] <- data.frame(
        row_id = row_id,
        engine = as.character(control_metrics$engine %||% "not_run"),
        n_solves = as.integer(control_metrics$n_solves %||% 0L),
        n_objective_updates = as.integer(
          control_metrics$n_objective_updates %||% 0L
        ),
        n_fallback = as.integer(control_metrics$n_fallback %||% 0L),
        reused_from_primary = isTRUE(all(reuse_mask)),
        n_reused_from_primary = as.integer(sum(reuse_mask)),
        reuse_fraction = mean(reuse_mask),
        shared_target_engine = isTRUE(ctrl$shared_target_engine),
        stringsAsFactors = FALSE
      )
    }
    rm(result)
    unlink(checkpoint_files[[i]], force = TRUE)
    invisible(gc(verbose = FALSE, full = FALSE))
  }
  if (anyDuplicated(observed_rows) || !setequal(observed_rows, row_ids)) {
    stop("Step 2 reaction batches did not score every target exactly once.",
         call. = FALSE)
  }
  unlink(payload_files, force = TRUE)
  rm(checkpoint_groups, checkpoint_files, tasks, batch_specs)
  invisible(gc(verbose = FALSE, full = TRUE))

  score <- rc_compass_score_from_penalty(penalty, feasible)
  score_rna_only <- if (!is.null(control_penalty)) {
    rc_compass_score_from_penalty(control_penalty, control_feasible)
  } else NULL
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
    penalty_rna_only = control_penalty,
    score_rna_only = score_rna_only,
    feasible_rna_only = control_feasible,
    evaluated_rna_only = control_evaluated,
    lp_diagnostics_rna_only = .rc_bind_frames_fill(control_diagnostics),
    step2_engine_metrics_rna_only =
      .rc_bind_frames_fill(control_step2_engine_metrics),
    rna_control_model_identical_reuse = control_reused,
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
      step2_dispatch = "model_scoped_reaction_batches",
      step2_parallel_tasks = as.integer(step2_task_count),
      step2_parallel_workers = as.integer(step2_workers),
      step2_worker_payload = paste(
        "file-backed compact S/lb/ub, model-specific penalties, target",
        "metadata and cached vmax; no Layer 1/full GEM/global closure export"
      ),
      step2_model_load_reuse = paste(
        "controller validates each union GEM once; every reaction batch loads",
        "only compact S/lb/ub plus its matching-unit penalty matrix"
      ),
      vmax_computation_scope = vmax_computation_scope,
      vmax_parallel_tasks = vmax_parallel_tasks,
      vmax_parallel_workers = vmax_parallel_workers,
      vmax_solve_count = vmax_solve_count,
      vmax_reuse_by_cell_type = stats::setNames(
        as.integer(reuse), names(reuse)
      ),
      step2_solver_reuse = paste(
        "one prepared target template with independent persistent HiGHS",
        "streams for primary and RNA-only objectives across matching metacells"
      ),
      paired_primary_rna_control = !is.null(control_penalty),
      rna_control_vmax_solve_count = if (!is.null(control_penalty)) 0L else NA_integer_,
      rna_control_identical_model_reuse = sum(control_reused),
      worker_cleanup = paste(
        "checkpoint each reaction; release its target HiGHS engine; retain",
        "only one compact S/lb/ub model and matching penalty payload per batch"
      ),
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
  progress_state <- get0(
    ".rc_layer2_progress_state", mode = "environment", inherits = TRUE
  )
  args <- list(...)
  answer <- do.call(
    .rc_run_celltype_microcompass_engine_reaction_core,
    args
  )
  contexts <- .rc_layer2_model_contexts(
    answer$shared_model_cache, mode = "meta_module_gem"
  )
  parts_dir <- .rc_layer2_cache_progress_dir(answer$shared_model_cache)
  run_kind <- progress_state$run_kind %||% "primary"
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
