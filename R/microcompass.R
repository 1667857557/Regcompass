# Two-mode microCOMPASS runner.

#' Prepare direction-specific reaction targets from signed Human-GEM bounds
#' @export
rc_prepare_directional_targets <- function(gem, target_reactions,
                                           target_direction = c(
                                             "both", "forward", "reverse"
                                           ),
                                           bound_tolerance = 1e-12) {
  target_direction <- match.arg(target_direction)
  if (!is.numeric(bound_tolerance) || length(bound_tolerance) != 1L ||
      !is.finite(bound_tolerance) || bound_tolerance < 0) {
    stop("`bound_tolerance` must be one finite non-negative number.", call. = FALSE)
  }
  validated <- rc_validate_gem(gem)
  requested <- unique(trimws(as.character(target_reactions)))
  requested <- requested[!is.na(requested) & nzchar(requested)]
  missing <- setdiff(requested, validated$reactions)
  if (length(missing)) {
    stop("Target reactions missing from GEM: ",
         paste(utils::head(missing, 10L), collapse = ", "), call. = FALSE)
  }
  if (!length(requested)) {
    return(data.frame(
      reaction_id = character(), target_direction = character(),
      direction_class = character(), requested_direction = character(),
      direction_status = character(), stringsAsFactors = FALSE
    ))
  }

  rows <- lapply(requested, function(reaction) {
    lb <- validated$lb[[reaction]]
    ub <- validated$ub[[reaction]]
    forward_allowed <- ub > bound_tolerance
    reverse_allowed <- lb < -bound_tolerance
    direction_class <- if (forward_allowed && reverse_allowed) {
      "reversible"
    } else if (forward_allowed) {
      "forward_only"
    } else if (reverse_allowed) {
      "reverse_only"
    } else {
      "blocked"
    }
    directions <- switch(
      target_direction,
      both = c(if (forward_allowed) "forward", if (reverse_allowed) "reverse"),
      forward = if (forward_allowed) "forward" else character(),
      reverse = if (reverse_allowed) "reverse" else character()
    )
    if (!length(directions)) directions <- "none"
    data.frame(
      reaction_id = reaction,
      target_direction = directions,
      direction_class = direction_class,
      requested_direction = target_direction,
      direction_status = ifelse(directions == "none", "no_allowed_direction", "allowed"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.rc_normalize_medium_scenarios <- function(medium_scenarios) {
  if (is.null(medium_scenarios) ||
      (is.data.frame(medium_scenarios) && !nrow(medium_scenarios))) {
    return(data.frame(
      medium_scenario_id = "base",
      exchange_reaction_id = NA_character_,
      lb = NA_real_,
      ub = NA_real_,
      available = FALSE,
      .no_constraints = TRUE,
      stringsAsFactors = FALSE
    ))
  }
  if (!is.data.frame(medium_scenarios)) {
    stop("`medium_scenarios` must be a data.frame.", call. = FALSE)
  }
  if (!"medium_scenario_id" %in% colnames(medium_scenarios)) {
    medium_scenarios$medium_scenario_id <- "custom"
  }
  medium_scenarios$.no_constraints <- FALSE
  medium_scenarios
}

.rc_cache_gem <- function(entry) {
  if (is.list(entry) && !is.null(entry$file)) {
    readRDS(entry$file)
  } else {
    entry
  }
}

.rc_unit_sample_map <- function(unit_meta, unit_ids, sample_col) {
  if (is.null(unit_meta) || !sample_col %in% colnames(unit_meta)) {
    stop(
      "Meta-module GEM mode requires unit metadata with `",
      sample_col,
      "`.",
      call. = FALSE
    )
  }
  id_columns <- intersect(
    c("unit_id", "pool_id", "metacell_id"),
    colnames(unit_meta)
  )
  for (id_column in id_columns) {
    output <- stats::setNames(
      as.character(unit_meta[[sample_col]]),
      as.character(unit_meta[[id_column]])
    )
    if (all(unit_ids %in% names(output))) {
      return(output[unit_ids])
    }
  }
  if (nrow(unit_meta) == length(unit_ids)) {
    return(stats::setNames(
      as.character(unit_meta[[sample_col]]),
      unit_ids
    ))
  }
  stop("Could not align Layer 2 units to sample IDs.", call. = FALSE)
}

.rc_unit_condition_map <- function(unit_meta, unit_ids, condition_col) {
  if (is.null(condition_col) || !nzchar(condition_col) ||
      is.null(unit_meta) || !condition_col %in% colnames(unit_meta)) {
    return(stats::setNames(rep("all", length(unit_ids)), unit_ids))
  }
  id_columns <- intersect(c("unit_id", "pool_id", "metacell_id"), colnames(unit_meta))
  for (id_column in id_columns) {
    values <- stats::setNames(
      as.character(unit_meta[[condition_col]]),
      as.character(unit_meta[[id_column]])
    )
    if (all(unit_ids %in% names(values))) {
      out <- values[unit_ids]
      if (anyNA(out) || any(!nzchar(out))) {
        stop("Unit conditions must be non-missing and non-empty.", call. = FALSE)
      }
      return(out)
    }
  }
  stop("Could not align Layer 2 units to conditions.", call. = FALSE)
}

.rc_sample_condition_map <- function(unit_meta, sample_col, condition_col) {
  if (is.null(unit_meta) || !sample_col %in% colnames(unit_meta)) {
    stop("Unit metadata is missing the sample column.", call. = FALSE)
  }
  samples <- unique(as.character(unit_meta[[sample_col]]))
  if (is.null(condition_col) || !nzchar(condition_col) ||
      !condition_col %in% colnames(unit_meta)) {
    return(stats::setNames(rep("all", length(samples)), samples))
  }
  pairs <- unique(data.frame(
    sample_id = as.character(unit_meta[[sample_col]]),
    condition = as.character(unit_meta[[condition_col]]),
    stringsAsFactors = FALSE
  ))
  if (anyNA(pairs$sample_id) || any(!nzchar(pairs$sample_id)) ||
      anyNA(pairs$condition) || any(!nzchar(pairs$condition)) ||
      any(table(pairs$sample_id) != 1L)) {
    stop("Each biological sample must map to exactly one condition.", call. = FALSE)
  }
  stats::setNames(pairs$condition, pairs$sample_id)
}


#' Run microCOMPASS with a full GEM or a FASTCORE-completed meta-module GEM
#'
#' The package intentionally exposes exactly two structural modes.
#' @export
rc_run_microcompass <- function(layer1, gem,
                                 target_reactions = NULL,
                                 medium_table = NULL,
                                 medium_scenarios = NULL,
                                 mode = c(
                                   "full_gem", "meta_module_gem"
                                 ),
                                 reaction_membership = NULL,
                                 core_reactions = NULL,
                                 unit = c(
                                   "sample_celltype", "metacell"
                                 ),
                                 condition_col = "condition",
                                 sample_col = "sample_id",
                                 celltype_col = "cell_type",
                                 model_params = list(),
                                 penalty_weights = c(
                                   expr = 1.0,
                                   confidence = 0.5,
                                   missing = 1.0
                                 ),
                                 omega = 0.95,
                                 target_direction = c(
                                   "both", "forward", "reverse"
                                 ),
                                 parallel = TRUE,
                                 solver = c(
                                   "highs", "gurobi", "glpk"
                                 ),
                                 time_limit = 60,
                                 flux_threshold = 1e-8,
                                 BPPARAM = NULL) {
  mode <- match.arg(mode)
  unit <- match.arg(unit)
  solver <- match.arg(solver)
  target_direction <- match.arg(target_direction)
  medium_scenarios <- .rc_normalize_medium_scenarios(
    medium_scenarios %||% medium_table
  )
  matrices <- rc_layer2_unit_matrices(
    layer1,
    if (unit == "metacell") "metacell" else "sample_celltype",
    sample_col,
    celltype_col,
    condition_col
  )
  gem <- rc_annotate_reaction_roles(gem)
  unit_condition <- .rc_unit_condition_map(
    matrices$unit_meta, colnames(matrices$C_rel), condition_col
  )
  sample_conditions <- .rc_sample_condition_map(
    matrices$unit_meta, sample_col, condition_col
  )
  direction_diagnostics <- NULL

  if (identical(mode, "full_gem")) {
    if (is.null(target_reactions) || !length(target_reactions)) {
      stop(
        "`target_reactions` is required in full-GEM mode.",
        call. = FALSE
      )
    }
    directions <- rc_prepare_directional_targets(
      gem,
      target_reactions,
      target_direction
    )
    direction_diagnostics <- directions
    directions <- directions[
      directions$target_direction %in% c("forward", "reverse"),
      , drop = FALSE
    ]
    if (!nrow(directions)) {
      stop(
        "No target reaction directions are allowed by the GEM bounds.",
        call. = FALSE
      )
    }
    cache_dir <- model_params$cache_dir %||%
      tempfile("RegCompassR_full_gem_cache_")
    model_cache <- rc_build_full_gem_cache(
      gem = gem,
      dirs = directions,
      medium_scenarios = medium_scenarios,
      cache_dir = cache_dir,
      conditions = unique(unname(unit_condition))
    )
  } else {
    if (is.null(reaction_membership)) {
      stop(
        paste(
          "`reaction_membership` is required in",
          "meta-module-GEM mode."
        ),
        call. = FALSE
      )
    }
    cache_dir <- model_params$cache_dir %||%
      tempfile("RegCompassR_meta_module_gem_cache_")
    merge_sample_modules <- isTRUE(
      model_params$merge_sample_modules %||% TRUE
    )
    cache_membership <- reaction_membership
    cache_core <- core_reactions
    cache_sample_conditions <- sample_conditions
    if (merge_sample_modules) {
      merged_sample_id <- model_params$merged_sample_id %||%
        "__merged_meta_module__"
      merged_module_id <- model_params$merged_module_id %||%
        "merged_meta_module"
      cache_membership$source_sample_id <- as.character(
        cache_membership$sample_id
      )
      cache_membership$source_module_id <- as.character(
        cache_membership$module_id
      )
      cache_membership$sample_id <- merged_sample_id
      cache_membership$module_id <- merged_module_id
      if (!is.null(cache_core)) {
        cache_core$source_sample_id <- as.character(cache_core$sample_id)
        cache_core$source_module_id <- as.character(cache_core$module_id)
        cache_core$sample_id <- merged_sample_id
        cache_core$module_id <- merged_module_id
      }
      cache_sample_conditions <- data.frame(
        sample_id = merged_sample_id,
        condition = "all",
        stringsAsFactors = FALSE
      )
    }
    model_cache <- rc_build_meta_module_gem_cache(
      gem = gem,
      reaction_membership = cache_membership,
      core_reactions = cache_core,
      target_reactions = target_reactions,
      medium_scenarios = medium_scenarios,
      cache_dir = cache_dir,
      target_direction = target_direction,
      solver = solver,
      time_limit = model_params$completion_time_limit %||%
        max(time_limit, 300),
      fastcore_epsilon = model_params$fastcore_epsilon %||%
        1e-4,
      max_support_reactions =
        model_params$max_support_reactions %||% 2000,
      strict = model_params$strict %||% TRUE,
      sample_conditions = cache_sample_conditions
    )
    attr(model_cache, "merge_sample_modules") <- merge_sample_modules
    if (!length(model_cache)) {
      stop(
        "No parent-feasible meta-module targets were available.",
        call. = FALSE
      )
    }
    directions <- unique(do.call(rbind, lapply(
      model_cache,
      function(entry) {
        data.frame(
          sample_id = entry$sample_id,
          module_id = entry$module_id,
          reaction_id = entry$reaction_id,
          target_direction = entry$target_direction,
          condition = entry$condition %||% "all",
          stringsAsFactors = FALSE
        )
      }
    )))
  }

  all_reactions <- unique(unlist(lapply(
    model_cache,
    function(entry) colnames(.rc_cache_gem(entry)$S)
  )))
  penalties <- rc_compute_multiome_penalty(
    rc_align_layer2_evidence(
      matrices$C_rel,
      all_reactions,
      NA_real_
    ),
    rc_align_layer2_evidence(
      matrices$reaction_confidence,
      all_reactions,
      NA_real_
    ),
    layer1$gpr_diagnostics,
    gem$reaction_roles,
    weights = penalty_weights
  )

  units <- colnames(matrices$C_rel)
  row_ids <- names(model_cache)
  penalty <- vmax <- matrix(
    NA_real_,
    nrow = length(row_ids),
    ncol = length(units),
    dimnames = list(row_ids, units)
  )
  feasible <- matrix(
    FALSE,
    nrow = length(row_ids),
    ncol = length(units),
    dimnames = list(row_ids, units)
  )
  evaluated <- feasible

  if (identical(mode, "meta_module_gem")) {
    merged_meta_module <- isTRUE(
      attr(model_cache, "merge_sample_modules")
    )
    unit_sample <- .rc_unit_sample_map(
      matrices$unit_meta,
      units,
      sample_col
    )
    task_rows <- lapply(row_ids, function(row_id) {
      entry <- model_cache[[row_id]]
      selected_units <- if (merged_meta_module) {
        condition <- entry$condition %||% "all"
        if (identical(condition, "all")) {
          units
        } else {
          names(unit_condition)[unit_condition == condition]
        }
      } else {
        names(unit_sample)[
          unit_sample == as.character(entry$sample_id)
        ]
      }
      if (!length(selected_units)) return(NULL)
      data.frame(
        row_id = row_id,
        unit_id = selected_units,
        stringsAsFactors = FALSE
      )
    })
    task_rows <- task_rows[
      !vapply(task_rows, is.null, logical(1))
    ]
    tasks <- if (length(task_rows)) {
      do.call(rbind, task_rows)
    } else {
      data.frame()
    }
  } else {
    task_rows <- lapply(row_ids, function(row_id) {
      entry <- model_cache[[row_id]]
      condition <- entry$condition %||% "all"
      selected_units <- if (identical(condition, "all")) {
        units
      } else {
        names(unit_condition)[unit_condition == condition]
      }
      if (!length(selected_units)) return(NULL)
      data.frame(
        row_id = row_id,
        unit_id = selected_units,
        stringsAsFactors = FALSE
      )
    })
    task_rows <- task_rows[!vapply(task_rows, is.null, logical(1))]
    tasks <- if (length(task_rows)) do.call(rbind, task_rows) else data.frame()
  }
  if (!nrow(tasks)) {
    stop(
      "No model-target and unit combinations were available.",
      call. = FALSE
    )
  }

  run_one <- function(task) {
    entry <- model_cache[[task$row_id]]
    model <- .rc_cache_gem(entry)
    unit_id <- as.character(task$unit_id)
    reaction_penalties <- penalties$penalty[
      colnames(model$S),
      unit_id
    ]
    answer <- rc_compass_two_step_lp_directional(
      S = model$S,
      lb = model$lb,
      ub = model$ub,
      target_reaction = entry$reaction_id,
      penalties = reaction_penalties,
      target_direction = entry$target_direction,
      omega = omega,
      solver = solver,
      time_limit = time_limit,
      flux_threshold = flux_threshold
    )
    diagnostic <- data.frame(
      row_id = task$row_id,
      unit_id = unit_id,
      sample_id = entry$sample_id %||% NA_character_,
      module_id = entry$module_id %||% NA_character_,
      reaction_id = entry$reaction_id,
      target_direction = entry$target_direction,
      medium_scenario = entry$medium_scenario,
      condition = entry$condition %||% "all",
      strict_feasible = isTRUE(answer$feasible),
      solver_status = answer$solver_status,
      step1_status = answer$step1_status,
      step2_status = answer$step2_status,
      target_status = model$target_status %||%
        if (isTRUE(answer$feasible)) {
          "ok"
        } else {
          "structurally_infeasible"
        },
      objective_value = answer$penalty,
      vmax = answer$vmax,
      stringsAsFactors = FALSE
    )
    list(
      row_id = task$row_id,
      unit_id = unit_id,
      penalty = answer$penalty,
      vmax = answer$vmax,
      feasible = isTRUE(answer$feasible),
      diagnostics = diagnostic
    )
  }

  task_list <- split(tasks, seq_len(nrow(tasks)))
  results <- rc_parallel_lapply(
    task_list,
    function(task) run_one(task[1, , drop = FALSE]),
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )
  for (result in results) {
    penalty[result$row_id, result$unit_id] <- result$penalty
    vmax[result$row_id, result$unit_id] <- result$vmax
    feasible[result$row_id, result$unit_id] <- result$feasible
    evaluated[result$row_id, result$unit_id] <- TRUE
  }
  score <- rc_compass_score_from_penalty(penalty, feasible)
  lp_diagnostics <- do.call(
    rbind,
    lapply(results, `[[`, "diagnostics")
  )

  diagnostic_keys <- vapply(
    model_cache,
    function(entry) {
      entry$file %||% paste(
        entry$sample_id %||% "",
        entry$module_id %||% "",
        entry$medium_scenario,
        sep = "::"
      )
    },
    character(1)
  )
  diagnostic_entries <- model_cache[!duplicated(diagnostic_keys)]
  diagnostic_rows <- lapply(diagnostic_entries, function(entry) {
    model <- .rc_cache_gem(entry)
    output <- model$closure_diagnostics
    if (!is.data.frame(output) || !nrow(output)) return(NULL)
    output$sample_id <- entry$sample_id %||% NA_character_
    output$module_id <- entry$module_id %||% NA_character_
    output$medium_scenario <- entry$medium_scenario
    output
  })
  diagnostic_rows <- diagnostic_rows[
    !vapply(diagnostic_rows, is.null, logical(1))
  ]
  model_diagnostics <- if (length(diagnostic_rows)) {
    do.call(rbind, diagnostic_rows)
  } else {
    data.frame()
  }

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
    unit_meta = matrices$unit_meta,
    params = list(
      unit = unit,
      omega = omega,
      target_direction = target_direction,
      model_mode = mode,
      merge_sample_modules = if (identical(mode, "meta_module_gem")) {
        isTRUE(attr(model_cache, "merge_sample_modules"))
      } else {
        FALSE
      },
      flux_threshold = flux_threshold,
      evidence_policy = paste(
        "RNA+ATAC-GPR evidence affects COMPASS penalties only;",
        "meta-module structure is biological membership plus",
        "FASTCORE support."
      )
    ),
    method = if (identical(mode, "full_gem")) {
      "microCOMPASS full-GEM directional LP"
    } else {
      paste(
        "microCOMPASS FASTCORE-completed",
        "meta-module-GEM directional LP"
      )
    }
  )
}

#' Summarize a microCOMPASS result
#' @export
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
    stringsAsFactors = FALSE
  )
}
