# Shared-structure directional microCOMPASS runner.

#' Prepare direction-specific reaction targets from active GEM bounds
#' @export
rc_prepare_directional_targets <- function(gem, target_reactions,
                                           target_direction = c("both", "forward", "reverse"),
                                           bound_tolerance = 1e-12) {
  target_direction <- match.arg(target_direction)
  validated <- rc_validate_gem(gem)
  requested <- unique(trimws(as.character(target_reactions)))
  requested <- requested[!is.na(requested) & nzchar(requested)]
  missing <- setdiff(requested, validated$reactions)
  if (length(missing)) stop("Target reactions missing from GEM: ", paste(utils::head(missing, 10L), collapse = ", "), call. = FALSE)
  rows <- lapply(requested, function(reaction) {
    forward <- validated$ub[[reaction]] > bound_tolerance
    reverse <- validated$lb[[reaction]] < -bound_tolerance
    allowed <- switch(target_direction,
      both = c(if (forward) "forward", if (reverse) "reverse"),
      forward = if (forward) "forward" else character(),
      reverse = if (reverse) "reverse" else character()
    )
    if (!length(allowed)) allowed <- "none"
    data.frame(
      reaction_id = reaction,
      target_direction = allowed,
      direction_class = if (forward && reverse) "reversible" else if (forward) "forward_only" else if (reverse) "reverse_only" else "blocked",
      requested_direction = target_direction,
      direction_status = ifelse(allowed == "none", "no_allowed_direction", "allowed"),
      stringsAsFactors = FALSE
    )
  })
  if (length(rows)) do.call(rbind, rows) else data.frame(
    reaction_id = character(), target_direction = character(), direction_class = character(),
    requested_direction = character(), direction_status = character(), stringsAsFactors = FALSE
  )
}

.rc_normalize_medium_scenarios <- function(medium_scenarios) {
  if (is.null(medium_scenarios) || (is.data.frame(medium_scenarios) && !nrow(medium_scenarios))) {
    return(data.frame(
      medium_scenario_id = "base", exchange_reaction_id = NA_character_, lb = NA_real_,
      ub = NA_real_, available = FALSE, .no_constraints = TRUE, stringsAsFactors = FALSE
    ))
  }
  if (!is.data.frame(medium_scenarios)) stop("`medium_scenarios` must be a data.frame.", call. = FALSE)
  if (!"medium_scenario_id" %in% colnames(medium_scenarios)) medium_scenarios$medium_scenario_id <- "custom"
  medium_scenarios$.no_constraints <- FALSE
  medium_scenarios
}

.rc_validate_shared_medium <- function(medium_scenarios) {
  if ("condition" %in% colnames(medium_scenarios)) {
    values <- unique(trimws(as.character(medium_scenarios$condition)))
    values <- values[!is.na(values) & nzchar(values)]
    if (length(setdiff(values, "all"))) {
      stop("Shared-GEM scoring requires condition-independent medium bounds; use condition = 'all'.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

.rc_cache_gem <- function(entry) if (is.list(entry) && !is.null(entry$file)) readRDS(entry$file) else entry

.rc_global_union_tables <- function(reaction_membership, core_reactions) {
  if (!is.data.frame(reaction_membership) || !"reaction_id" %in% colnames(reaction_membership)) {
    stop("`reaction_membership` must contain `reaction_id`.", call. = FALSE)
  }
  if (is.null(core_reactions)) core_reactions <- reaction_membership
  if (!is.data.frame(core_reactions) || !"reaction_id" %in% colnames(core_reactions)) {
    stop("`core_reactions` must contain `reaction_id`.", call. = FALSE)
  }
  membership <- unique(data.frame(
    sample_id = "GLOBAL", module_id = "GLOBAL_UNION",
    reaction_id = trimws(as.character(reaction_membership$reaction_id)),
    is_core = FALSE, stringsAsFactors = FALSE
  ))
  membership <- membership[!is.na(membership$reaction_id) & nzchar(membership$reaction_id), , drop = FALSE]
  hard <- if ("is_core" %in% colnames(core_reactions)) core_reactions$is_core %in% TRUE else rep(TRUE, nrow(core_reactions))
  core_ids <- unique(trimws(as.character(core_reactions$reaction_id[hard])))
  core_ids <- core_ids[!is.na(core_ids) & nzchar(core_ids)]
  membership$is_core <- membership$reaction_id %in% core_ids
  core <- membership[membership$is_core, , drop = FALSE]
  if (!nrow(core)) stop("No hard core reactions remain after global union.", call. = FALSE)
  list(membership = membership, core = core)
}

.rc_build_global_meta_module_cache <- function(gem, reaction_membership, core_reactions,
                                               target_reactions, medium_scenarios,
                                               cache_dir, target_direction, solver,
                                               time_limit, model_params) {
  union <- .rc_global_union_tables(reaction_membership, core_reactions)
  medium_scenarios <- .rc_normalize_medium_scenarios(medium_scenarios)
  .rc_validate_shared_medium(medium_scenarios)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  scenarios <- unique(as.character(medium_scenarios$medium_scenario_id))
  cache <- list(); summaries <- list(); direction_rows <- list(); index <- 1L
  for (scenario in scenarios) {
    medium <- medium_scenarios[as.character(medium_scenarios$medium_scenario_id) == scenario, , drop = FALSE]
    if (!nrow(medium) || (".no_constraints" %in% colnames(medium) && all(medium$.no_constraints))) medium <- NULL
    file <- file.path(cache_dir, paste0("global_union__medium_", gsub("[^A-Za-z0-9_.-]+", "_", scenario), ".rds"))
    model <- rc_build_meta_module_gem(
      gem = gem, reaction_membership = union$membership, core_reactions = union$core,
      sample_id = "GLOBAL", module_id = "GLOBAL_UNION", medium_table = medium,
      condition = NULL, target_direction = target_direction, solver = solver,
      time_limit = model_params$completion_time_limit %||% max(time_limit, 300),
      fastcore_epsilon = model_params$fastcore_epsilon %||% 1e-4,
      max_support_reactions = model_params$max_support_reactions %||% 2000,
      strict = model_params$strict %||% TRUE
    )
    model$global_union <- TRUE
    saveRDS(model, file)
    allowed <- model$target_directions
    if (!is.null(target_reactions)) allowed <- allowed[allowed$reaction_id %in% target_reactions, , drop = FALSE]
    if (nrow(allowed)) {
      for (j in seq_len(nrow(allowed))) {
        key <- paste0(
          "reaction=", utils::URLencode(as.character(allowed$reaction_id[[j]]), reserved = TRUE),
          "::direction=", as.character(allowed$target_direction[[j]]),
          "::medium=", utils::URLencode(scenario, reserved = TRUE), "::condition=all"
        )
        cache[[key]] <- list(
          reaction_id = as.character(allowed$reaction_id[[j]]),
          target_direction = as.character(allowed$target_direction[[j]]),
          medium_scenario = scenario, condition = "all", file = file,
          build_strategy = "global_union_meta_module_gem"
        )
      }
      allowed$medium_scenario <- scenario
      direction_rows[[length(direction_rows) + 1L]] <- allowed
    }
    summaries[[index]] <- data.frame(
      cache_key = paste("global_union", scenario, sep = "::"), medium_scenario = scenario,
      condition = "all", file = file, n_reactions = ncol(model$S), n_metabolites = nrow(model$S),
      n_biological_reactions = model$build_params$n_biological_reactions,
      n_fastcore_support_reactions = model$build_params$n_fastcore_support_reactions,
      target_status = model$target_status, build_strategy = "global_union_meta_module_gem",
      stringsAsFactors = FALSE
    )
    index <- index + 1L
  }
  attr(cache, "summary") <- do.call(rbind, summaries)
  attr(cache, "target_directions") <- if (length(direction_rows)) do.call(rbind, direction_rows) else data.frame()
  attr(cache, "global_union") <- union
  cache
}

#' Run directional microCOMPASS with one shared structural model per medium
#' @export
rc_run_microcompass <- function(layer1, gem, target_reactions = NULL,
                                 medium_table = NULL, medium_scenarios = NULL,
                                 mode = c("full_gem", "meta_module_gem"),
                                 reaction_membership = NULL, core_reactions = NULL,
                                 unit = c("metacell", "sample_celltype"),
                                 condition_col = "condition", sample_col = "sample_id",
                                 celltype_col = "cell_type", model_params = list(),
                                 penalty_weights = c(expr = 1, confidence = 0.5, missing = 1),
                                 omega = 0.95,
                                 target_direction = c("both", "forward", "reverse"),
                                 parallel = TRUE, solver = c("highs", "gurobi", "glpk"),
                                 time_limit = 60, flux_threshold = 1e-8, BPPARAM = NULL) {
  mode <- match.arg(mode); unit <- match.arg(unit); solver <- match.arg(solver); target_direction <- match.arg(target_direction)
  medium_scenarios <- .rc_normalize_medium_scenarios(medium_scenarios %||% medium_table)
  .rc_validate_shared_medium(medium_scenarios)
  matrices <- rc_layer2_unit_matrices(layer1, unit, sample_col, celltype_col, condition_col)
  gem <- rc_annotate_reaction_roles(gem)
  cache_dir <- model_params$cache_dir %||% tempfile("RegCompassR_shared_gem_cache_")
  if (identical(mode, "full_gem")) {
    if (is.null(target_reactions) || !length(target_reactions)) stop("`target_reactions` is required in full-GEM mode.", call. = FALSE)
    diagnostics <- rc_prepare_directional_targets(gem, target_reactions, target_direction)
    directions <- diagnostics[diagnostics$target_direction %in% c("forward", "reverse"), , drop = FALSE]
    if (!nrow(directions)) stop("No target reaction directions are allowed by the GEM bounds.", call. = FALSE)
    model_cache <- rc_build_full_gem_cache(gem, directions, medium_scenarios, cache_dir = cache_dir, conditions = "all")
  } else {
    if (is.null(reaction_membership)) stop("`reaction_membership` is required in meta-module-GEM mode.", call. = FALSE)
    model_cache <- .rc_build_global_meta_module_cache(
      gem, reaction_membership, core_reactions, target_reactions, medium_scenarios,
      cache_dir, target_direction, solver, time_limit, model_params
    )
    directions <- attr(model_cache, "target_directions")
    diagnostics <- directions
  }
  if (!length(model_cache)) stop("No shared-GEM targets were available for scoring.", call. = FALSE)
  all_reactions <- unique(unlist(lapply(model_cache, function(entry) colnames(.rc_cache_gem(entry)$S))))
  penalties <- rc_compute_multiome_penalty(
    rc_align_layer2_evidence(matrices$C_rel, all_reactions, NA_real_),
    rc_align_layer2_evidence(matrices$reaction_confidence, all_reactions, NA_real_),
    layer1$gpr_diagnostics, gem$reaction_roles, weights = penalty_weights
  )
  units <- colnames(matrices$C_rel); row_ids <- names(model_cache)
  penalty <- vmax <- matrix(NA_real_, length(row_ids), length(units), dimnames = list(row_ids, units))
  feasible <- evaluated <- matrix(FALSE, length(row_ids), length(units), dimnames = list(row_ids, units))
  tasks <- expand.grid(row_id = row_ids, unit_id = units, stringsAsFactors = FALSE)
  run_one <- function(task) {
    entry <- model_cache[[as.character(task$row_id)]]
    model <- .rc_cache_gem(entry)
    unit_id <- as.character(task$unit_id)
    answer <- rc_compass_two_step_lp_directional(
      S = model$S, lb = model$lb, ub = model$ub,
      target_reaction = entry$reaction_id,
      penalties = penalties$penalty[colnames(model$S), unit_id],
      target_direction = entry$target_direction, omega = omega,
      solver = solver, time_limit = time_limit, flux_threshold = flux_threshold
    )
    list(
      row_id = as.character(task$row_id), unit_id = unit_id,
      penalty = answer$penalty, vmax = answer$vmax, feasible = isTRUE(answer$feasible),
      diagnostics = data.frame(
        row_id = as.character(task$row_id), unit_id = unit_id,
        reaction_id = entry$reaction_id, target_direction = entry$target_direction,
        medium_scenario = entry$medium_scenario, condition = "all",
        strict_feasible = isTRUE(answer$feasible), solver_status = answer$solver_status,
        step1_status = answer$step1_status, step2_status = answer$step2_status,
        objective_value = answer$penalty, vmax = answer$vmax, stringsAsFactors = FALSE
      )
    )
  }
  task_list <- split(tasks, seq_len(nrow(tasks)))
  results <- rc_parallel_lapply(task_list, function(x) run_one(x[1, , drop = FALSE]), BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE)
  for (answer in results) {
    penalty[answer$row_id, answer$unit_id] <- answer$penalty
    vmax[answer$row_id, answer$unit_id] <- answer$vmax
    feasible[answer$row_id, answer$unit_id] <- answer$feasible
    evaluated[answer$row_id, answer$unit_id] <- TRUE
  }
  score <- rc_compass_score_from_penalty(penalty, feasible)
  model_files <- unique(vapply(model_cache, function(entry) as.character(entry$file), character(1)))
  diagnostic_rows <- lapply(model_files, function(file) {
    model <- readRDS(file)
    diagnostics <- model$closure_diagnostics
    if (!is.data.frame(diagnostics) || !nrow(diagnostics)) return(NULL)
    diagnostics$model_file <- file
    diagnostics
  })
  diagnostic_rows <- diagnostic_rows[!vapply(diagnostic_rows, is.null, logical(1))]
  model_diagnostics <- if (length(diagnostic_rows)) do.call(rbind, diagnostic_rows) else data.frame()
  list(
    score = score, penalty = penalty, vmax = vmax, feasible = feasible, evaluated = evaluated,
    target_direction = directions, direction_diagnostics = diagnostics,
    medium_scenarios = medium_scenarios, model_mode = mode,
    model_cache_summary = attr(model_cache, "summary"),
    model_diagnostics = model_diagnostics,
    lp_diagnostics = do.call(rbind, lapply(results, `[[`, "diagnostics")),
    penalty_components = penalties$components, evidence_policy = penalties$evidence_policy,
    unit_meta = matrices$unit_meta,
    params = list(unit = unit, omega = omega, target_direction = target_direction,
                  model_mode = mode, shared_structural_model = TRUE,
                  flux_threshold = flux_threshold),
    method = if (identical(mode, "full_gem")) "microCOMPASS shared full-GEM directional LP" else "microCOMPASS shared global-union meta-module GEM directional LP"
  )
}

#' Summarize a microCOMPASS result
#' @export
rc_summarize_microcompass <- function(result) {
  evaluated <- result$evaluated %||% matrix(TRUE, nrow(result$feasible), ncol(result$feasible), dimnames = dimnames(result$feasible))
  data.frame(
    model_mode = result$model_mode %||% NA_character_, n_targets = nrow(result$score),
    n_units = ncol(result$score), n_evaluated = sum(evaluated),
    feasible_fraction = if (any(evaluated)) mean(result$feasible[evaluated]) else NA_real_,
    stringsAsFactors = FALSE
  )
}
