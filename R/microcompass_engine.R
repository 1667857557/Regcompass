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

.rc_run_shared_full_gem_engine <- function(
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
  } else if (identical(mode, "full_gem")) {
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
  } else {
    stop("The shared full-GEM engine received a non-full-GEM mode.",
         call. = FALSE)
  }

  units <- colnames(matrices$reaction_expression)
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
  tasks <- expand.grid(
    model_key = unique_model_keys,
    unit_id = units,
    stringsAsFactors = FALSE
  )

  run_one_unit <- function(task) {
    model_key <- as.character(task$model_key)
    unit_id <- as.character(task$unit_id)
    selected_rows <- names(model_keys)[model_keys == model_key]
    first_entry <- model_cache[[selected_rows[[1L]]]]
    model <- .rc_load_microcompass_model(first_entry, mode)
    target_results <- lapply(selected_rows, function(row_id) {
      entry <- model_cache[[row_id]]
      unit_penalty <- penalties$penalty[colnames(model$S), unit_id]
      target_index <- match(entry$reaction_id, colnames(model$S))
      if (is.na(target_index)) {
        stop("A target reaction is absent from its shared model.",
             call. = FALSE)
      }
      target_penalty <- unit_penalty[[target_index]]
      evidence_available <- is.finite(unit_penalty)
      solver_penalty <- unit_penalty
      solver_penalty[!evidence_available] <- 0
      answer <- .rc_compass_step2_from_vmax_directional(
        S = model$S,
        lb = model$lb,
        ub = model$ub,
        target_reaction = entry$reaction_id,
        penalties = solver_penalty,
        vmax_result = vmax_cache[[row_id]],
        target_direction = entry$target_direction,
        omega = omega,
        solver = solver,
        flux_threshold = flux_threshold
      )
      target_evidence_available <- is.finite(target_penalty)
      list(
        row_id = row_id,
        unit_id = unit_id,
        penalty = if (target_evidence_available) {
          answer$penalty
        } else {
          NA_real_
        },
        vmax = answer$vmax,
        feasible = isTRUE(answer$feasible),
        evaluated = isTRUE(answer$feasible) &&
          target_evidence_available,
        diagnostics = data.frame(
          row_id = row_id,
          unit_id = unit_id,
          module_id = NA_character_,
          reaction_id = entry$reaction_id,
          target_direction = entry$target_direction,
          medium_scenario = entry$medium_scenario,
          condition = "all",
          strict_feasible = isTRUE(answer$feasible),
          solver_status = answer$solver_status,
          step1_status = answer$step1_status,
          step2_status = answer$step2_status,
          target_status = if (isTRUE(answer$feasible)) {
            "ok"
          } else {
            "medium_directionally_infeasible"
          },
          objective_value = if (target_evidence_available) {
            answer$penalty
          } else {
            NA_real_
          },
          vmax = answer$vmax,
          vmax_reused_from_shared_cache = TRUE,
          target_expression_available = target_evidence_available,
          objective_evidence_fraction = mean(evidence_available),
          unavailable_objective_terms = sum(!evidence_available),
          stringsAsFactors = FALSE
        )
      )
    })
    list(
      results = target_results,
      diagnostics = do.call(
        rbind, lapply(target_results, `[[`, "diagnostics")
      )
    )
  }

  task_list <- split(tasks, seq_len(nrow(tasks)))
  grouped_results <- rc_parallel_lapply(
    task_list,
    function(task) run_one_unit(task[1, , drop = FALSE]),
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )
  results <- unlist(
    lapply(grouped_results, `[[`, "results"), recursive = FALSE
  )
  for (result in results) {
    penalty[result$row_id, result$unit_id] <- result$penalty
    vmax[result$row_id, result$unit_id] <- result$vmax
    feasible[result$row_id, result$unit_id] <- result$feasible
    evaluated[result$row_id, result$unit_id] <- result$evaluated
  }
  score <- rc_compass_score_from_penalty(penalty, feasible)
  lp_diagnostics <- do.call(
    rbind, lapply(grouped_results, `[[`, "diagnostics")
  )
  model_diagnostics <- .rc_bind_frames_fill(lapply(
    representative_rows,
    function(row_id) {
      entry <- model_cache[[row_id]]
      model <- .rc_load_microcompass_model(entry, mode)
      model$closure_diagnostics %||% data.frame()
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
    penalty_components = penalties$components,
    evidence_policy = penalties$evidence_policy,
    evidence_policy_detail = penalties$evidence_policy_detail,
    unit_meta = matrices$unit_meta,
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
      parallel_task = "shared_model_by_metacell_step2",
      vmax_computation_scope =
        "shared_model_x_directional_target_once",
      vmax_solve_count = length(vmax_cache),
      vmax_reuse_factor = length(units),
      flux_threshold = flux_threshold,
      scoring_time_limit = "none"
    ),
    method = "microCOMPASS shared medium-constrained full-GEM directional LP"
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
