# Directional microCOMPASS scoring on shared structural models.

.rc_load_microcompass_model <- function(entry, mode) {
  if (identical(mode, "meta_module_gem")) {
    return(.rc_read_cached_union_gem(
      file = entry$file,
      medium_scenario = entry$medium_scenario,
      expected_checksum = entry$file_checksum %||% NA_character_
    ))
  }
  .rc_cache_gem(entry)
}

.rc_run_microcompass_engine <- function(
    layer1, gem, target_reactions = NULL,
    medium_table = NULL, medium_scenarios = NULL,
    mode = c("full_gem", "meta_module_gem"),
    reaction_membership = NULL, core_reactions = NULL,
    unit = c("metacell", "sample_celltype"),
    condition_col = "condition", sample_col = "sample_id",
    celltype_col = "cell_type", model_params = list(),
    omega = 0.95,
    target_direction = c("both", "forward", "reverse"),
    parallel = TRUE,
    solver = c("highs", "gurobi", "glpk"),
    flux_threshold = 1e-8,
    BPPARAM = NULL) {
  mode <- match.arg(mode)
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

  if (identical(mode, "full_gem")) {
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
      conditions = "all"
    )
  } else {
    model_cache <- .rc_build_medium_specific_union_gem_cache(
      gem = gem,
      reaction_membership = reaction_membership,
      core_reactions = core_reactions,
      target_reactions = target_reactions,
      medium_scenarios = medium_scenarios,
      cache_dir = model_params$cache_dir %||%
        tempfile("RegCompassR_medium_union_gem_cache_"),
      target_direction = target_direction,
      solver = solver,
      time_limit = model_params$completion_time_limit %||% 300,
      fastcore_epsilon = model_params$fastcore_epsilon %||% 1e-4,
      max_support_reactions = model_params$max_support_reactions %||% 2000,
      strict = model_params$strict %||% TRUE
    )
    if (!length(model_cache)) {
      stop("No parent-feasible union-GEM targets were available.",
           call. = FALSE)
    }
    directions <- unique(do.call(rbind, lapply(model_cache, function(entry) {
      data.frame(
        reaction_id = entry$reaction_id,
        target_direction = entry$target_direction,
        medium_scenario = entry$medium_scenario,
        stringsAsFactors = FALSE
      )
    })))
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
      answer <- rc_compass_two_step_lp_directional(
        S = model$S,
        lb = model$lb,
        ub = model$ub,
        target_reaction = entry$reaction_id,
        penalties = penalties$penalty[colnames(model$S), unit_id],
        target_direction = entry$target_direction,
        omega = omega,
        solver = solver,
        flux_threshold = flux_threshold
      )
      list(
        row_id = row_id,
        unit_id = unit_id,
        penalty = answer$penalty,
        vmax = answer$vmax,
        feasible = isTRUE(answer$feasible),
        diagnostics = data.frame(
          row_id = row_id,
          unit_id = unit_id,
          module_id = if (identical(mode, "meta_module_gem")) {
            "MEDIUM_UNION_GEM"
          } else {
            NA_character_
          },
          reaction_id = entry$reaction_id,
          target_direction = entry$target_direction,
          medium_scenario = entry$medium_scenario,
          condition = "all",
          strict_feasible = isTRUE(answer$feasible),
          solver_status = answer$solver_status,
          step1_status = answer$step1_status,
          step2_status = answer$step2_status,
          target_status = model$target_status %||%
            if (isTRUE(answer$feasible)) "ok" else "structurally_infeasible",
          objective_value = answer$penalty,
          vmax = answer$vmax,
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
    evaluated[result$row_id, result$unit_id] <- TRUE
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
    model_cache_summary = attr(model_cache, "summary"),
    model_diagnostics = model_diagnostics,
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
      shared_gem_scope = if (identical(mode, "meta_module_gem")) {
        "one_final_union_gem_per_medium_shared_across_all_units"
      } else {
        "one_full_gem_per_medium_shared_across_all_units"
      },
      parallel_task = "shared_model_by_metacell",
      flux_threshold = flux_threshold,
      scoring_time_limit = "none"
    ),
    method = if (identical(mode, "full_gem")) {
      "microCOMPASS shared full-GEM directional LP"
    } else {
      paste(
        "microCOMPASS directional LP on final medium-specific union GEMs",
        "after one global FASTCORE completion"
      )
    }
  )
}

#' Run directional minimum-evidence-discordance LPs
rc_run_microcompass <- function(
    layer1, gem, target_reactions = NULL,
    medium_table = NULL, medium_scenarios = NULL,
    mode = c("full_gem", "meta_module_gem"),
    reaction_membership = NULL, core_reactions = NULL,
    unit = c("sample_celltype", "metacell"),
    condition_col = "condition", sample_col = "sample_id",
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
        "Metacell-level scores are descriptive pseudo-observations and are not",
        "independent biological replicates."
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
      medium_scenarios <- rc_make_medium_scenarios(
        gem,
        scenario = "physiologic",
        species = "auto",
        strict_preset_matching = TRUE
      )
    } else {
      warning(
        paste(
          "GEM species provenance is unavailable; preserving the model's",
          "original exchange directions instead of assuming a human",
          "physiological medium."
        ),
        call. = FALSE
      )
      medium_scenarios <- rc_make_medium_scenarios(
        gem, scenario = "compass_model_bounds", species = "auto"
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
rc_summarize_microcompass <- function(result) {
  evaluated <- result$evaluated %||% matrix(
    TRUE,
    nrow = nrow(result$feasible),
    ncol = ncol(result$feasible),
    dimnames = dimnames(result$feasible)
  )
  data.frame(
    model_mode = result$model_mode %||% NA_character_,
    n_targets = nrow(result$score),
    n_units = ncol(result$score),
    n_evaluated = sum(evaluated),
    feasible_fraction = if (any(evaluated)) {
      mean(result$feasible[evaluated])
    } else {
      NA_real_
    },
    shared_gem = isTRUE(result$params$shared_gem),
    stringsAsFactors = FALSE
  )
}
