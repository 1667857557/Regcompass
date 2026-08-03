# Directional LP scoring on cell-type-specific structural models.

.rc_celltype_model_contract <- function(model_cache) {
  if (!is.list(model_cache) || !length(model_cache)) {
    stop("The cell-type model cache is empty.", call. = FALSE)
  }
  files <- vapply(model_cache, function(entry) as.character(entry$file), character(1))
  representative <- match(unique(files), files)
  rows <- lapply(representative, function(i) {
    entry <- model_cache[[i]]
    model <- .rc_read_celltype_union_gem(
      entry$file, entry$cell_type, entry$medium_scenario,
      entry$file_checksum
    )
    validated <- rc_validate_gem(model)
    data.frame(
      cell_type = as.character(entry$cell_type),
      medium_scenario = as.character(entry$medium_scenario),
      model_file = as.character(entry$file),
      model_file_checksum = unname(tools::md5sum(entry$file)[[1L]]),
      n_reactions = ncol(validated$S),
      n_metabolites = nrow(validated$S),
      reaction_order_checksum =
        .rc_microcompass_object_checksum(colnames(validated$S)),
      metabolite_order_checksum =
        .rc_microcompass_object_checksum(rownames(validated$S)),
      stoichiometry_bounds_checksum =
        .rc_microcompass_object_checksum(list(
          S = validated$S, lb = validated$lb, ub = validated$ub
        )),
      shared_across_conditions = TRUE,
      shared_across_cell_types = FALSE,
      stringsAsFactors = FALSE
    )
  })
  answer <- do.call(rbind, rows)
  rownames(answer) <- NULL
  answer[order(answer$cell_type, answer$medium_scenario), , drop = FALSE]
}

.rc_run_celltype_microcompass_engine <- function(
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
  unit_meta <- unit_meta[match(colnames(matrices$reaction_expression), unit_id),
                         , drop = FALSE]
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
  if (!setequal(cache_celltypes, unit_celltypes)) {
    stop(
      "Cell-type union GEMs and Layer 1 units cover different cell types.",
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
        "celltype_model_x_directional_target_once",
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
  task_rows <- lapply(unique_model_keys, function(model_key) {
    row_id <- representative_rows[[model_key]]
    cell_type <- as.character(model_cache[[row_id]]$cell_type)
    eligible <- names(unit_celltype)[unit_celltype == cell_type]
    if (!length(eligible)) {
      stop("No Layer 1 units match cell type `", cell_type, "`.",
           call. = FALSE)
    }
    data.frame(
      model_key = model_key,
      unit_id = eligible,
      stringsAsFactors = FALSE
    )
  })
  tasks <- do.call(rbind, task_rows)

  run_one_unit <- function(task) {
    model_key <- as.character(task$model_key)
    unit_id <- as.character(task$unit_id)
    selected_rows <- names(model_keys)[model_keys == model_key]
    first_entry <- model_cache[[selected_rows[[1L]]]]
    if (!identical(unit_celltype[[unit_id]], first_entry$cell_type)) {
      stop("A union GEM was assigned to a different cell type.", call. = FALSE)
    }
    model <- .rc_read_celltype_union_gem(
      first_entry$file, first_entry$cell_type,
      first_entry$medium_scenario, first_entry$file_checksum
    )
    target_results <- lapply(selected_rows, function(row_id) {
      entry <- model_cache[[row_id]]
      unit_penalty <- penalties$penalty[colnames(model$S), unit_id]
      target_index <- match(entry$reaction_id, colnames(model$S))
      if (is.na(target_index)) {
        stop("A target reaction is absent from its cell-type union GEM.",
             call. = FALSE)
      }
      target_penalty <- unit_penalty[[target_index]]
      evidence_available <- is.finite(unit_penalty)
      solver_penalty <- unit_penalty
      solver_penalty[!evidence_available] <- 0
      fit <- .rc_compass_step2_from_vmax_directional(
        S = model$S, lb = model$lb, ub = model$ub,
        target_reaction = entry$reaction_id,
        penalties = solver_penalty,
        vmax_result = vmax_cache[[row_id]],
        target_direction = entry$target_direction,
        omega = omega, solver = solver,
        flux_threshold = flux_threshold
      )
      target_available <- is.finite(target_penalty)
      list(
        row_id = row_id,
        unit_id = unit_id,
        penalty = if (target_available) fit$penalty else NA_real_,
        vmax = fit$vmax,
        feasible = isTRUE(fit$feasible),
        evaluated = isTRUE(fit$feasible) && target_available,
        diagnostics = data.frame(
          row_id = row_id,
          unit_id = unit_id,
          module_id = "CELLTYPE_MEDIUM_UNION_GEM",
          cell_type = entry$cell_type,
          reaction_id = entry$reaction_id,
          target_direction = entry$target_direction,
          medium_scenario = entry$medium_scenario,
          condition = "all",
          strict_feasible = isTRUE(fit$feasible),
          solver_status = fit$solver_status,
          step1_status = fit$step1_status,
          step2_status = fit$step2_status,
          target_status = model$target_status %||%
            if (isTRUE(fit$feasible)) "ok" else "structurally_infeasible",
          objective_value = if (target_available) fit$penalty else NA_real_,
          vmax = fit$vmax,
          vmax_reused_from_celltype_cache = TRUE,
          target_expression_available = target_available,
          objective_evidence_fraction = mean(evidence_available),
          unavailable_objective_terms = sum(!evidence_available),
          stringsAsFactors = FALSE
        )
      )
    })
    list(
      results = target_results,
      diagnostics = do.call(rbind, lapply(target_results, `[[`, "diagnostics"))
    )
  }

  grouped <- rc_parallel_lapply(
    split(tasks, seq_len(nrow(tasks))),
    function(task) run_one_unit(task[1, , drop = FALSE]),
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )
  results <- unlist(lapply(grouped, `[[`, "results"), recursive = FALSE)
  for (result in results) {
    penalty[result$row_id, result$unit_id] <- result$penalty
    vmax[result$row_id, result$unit_id] <- result$vmax
    feasible[result$row_id, result$unit_id] <- result$feasible
    evaluated[result$row_id, result$unit_id] <- result$evaluated
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
    lp_diagnostics = .rc_bind_frames_fill(lapply(grouped, `[[`, "diagnostics")),
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
      parallel_task = "celltype_model_by_matching_metacell_step2",
      vmax_computation_scope =
        "celltype_model_x_directional_target_once",
      vmax_solve_count = length(vmax_cache),
      vmax_reuse_by_cell_type = stats::setNames(
        as.integer(reuse), names(reuse)
      ),
      flux_threshold = flux_threshold,
      scoring_time_limit = "none"
    ),
    method = paste(
      "microCOMPASS directional LP on cell-type-specific medium union GEMs",
      "after independent FASTCORE completion within each cell type"
    )
  )
}
