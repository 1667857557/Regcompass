#' Prepare direction-specific reaction targets from signed GEM bounds
#' @export
rc_prepare_directional_targets <- function(gem, target_reactions,
                                           target_direction = c("both", "forward", "reverse"),
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
  if (length(missing)) stop("Target reactions missing from GEM: ", paste(utils::head(missing, 10L), collapse = ", "), call. = FALSE)
  rows <- lapply(requested, function(reaction) {
    lb <- validated$lb[[reaction]]
    ub <- validated$ub[[reaction]]
    forward_allowed <- ub > bound_tolerance
    reverse_allowed <- lb < -bound_tolerance
    direction_class <- if (forward_allowed && reverse_allowed) "reversible" else if (forward_allowed) "forward_only" else if (reverse_allowed) "reverse_only" else "blocked"
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
  if (!length(rows)) {
    return(data.frame(reaction_id = character(), target_direction = character(),
                      direction_class = character(), requested_direction = character(),
                      direction_status = character(), stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

.rc_normalize_medium_scenarios <- function(medium_scenarios) {
  if (is.null(medium_scenarios) || (is.data.frame(medium_scenarios) && !nrow(medium_scenarios))) {
    return(data.frame(
      medium_scenario_id = "base",
      exchange_reaction_id = NA_character_,
      lb = NA_real_, ub = NA_real_, available = FALSE,
      .no_constraints = TRUE, stringsAsFactors = FALSE
    ))
  }
  if (!is.data.frame(medium_scenarios)) stop("`medium_scenarios` must be a data.frame.", call. = FALSE)
  if (!"medium_scenario_id" %in% colnames(medium_scenarios)) medium_scenarios$medium_scenario_id <- "custom"
  medium_scenarios$.no_constraints <- FALSE
  medium_scenarios
}

.rc_validate_shared_medium <- function(medium_scenarios) {
  medium_scenarios <- .rc_normalize_medium_scenarios(medium_scenarios)
  if ("condition" %in% colnames(medium_scenarios)) {
    condition <- trimws(as.character(medium_scenarios$condition))
    condition <- unique(condition[!is.na(condition) & nzchar(condition) & condition != "all"])
    if (length(condition)) {
      stop("Shared-GEM scoring requires condition-invariant medium constraints; remove condition-specific medium rows.", call. = FALSE)
    }
  }
  medium_scenarios
}

.rc_cache_gem <- function(entry) {
  if (is.list(entry) && !is.null(entry$file)) readRDS(entry$file) else entry
}

.rc_build_global_meta_module_gem_cache <- function(gem, reaction_membership, core_reactions,
                                                    target_reactions = NULL, medium_scenarios = NULL,
                                                    cache_dir = tempfile("RegCompassR_global_meta_module_"),
                                                    target_direction = c("both", "forward", "reverse"),
                                                    solver = "highs", time_limit = 300,
                                                    fastcore_epsilon = 1e-4,
                                                    max_support_reactions = 2000,
                                                    strict = TRUE) {
  target_direction <- match.arg(target_direction)
  required <- c("sample_id", "module_id", "reaction_id")
  if (!is.data.frame(reaction_membership) || !all(required %in% colnames(reaction_membership))) {
    stop("`reaction_membership` must contain sample_id, module_id and reaction_id.", call. = FALSE)
  }
  if (!is.data.frame(core_reactions) || !all(required %in% colnames(core_reactions))) {
    stop("`core_reactions` must contain sample_id, module_id and reaction_id.", call. = FALSE)
  }
  if ("is_core" %in% colnames(core_reactions)) core_reactions <- core_reactions[core_reactions$is_core %in% TRUE, , drop = FALSE]
  if (!is.null(target_reactions)) {
    core_reactions <- core_reactions[as.character(core_reactions$reaction_id) %in% as.character(target_reactions), , drop = FALSE]
  }
  if (!nrow(core_reactions)) stop("No global core reactions remain for scoring.", call. = FALSE)
  reaction_membership$sample_id <- "global"
  reaction_membership$module_id <- "GLOBAL_UNION"
  core_reactions$sample_id <- "global"
  core_reactions$module_id <- "GLOBAL_UNION"
  medium_scenarios <- .rc_validate_shared_medium(medium_scenarios)
  scenarios <- unique(as.character(medium_scenarios$medium_scenario_id))
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  safe <- function(value) paste(sprintf("%02x", as.integer(charToRaw(enc2utf8(value)))), collapse = "")
  cache <- list()
  summaries <- list()
  for (scenario in scenarios) {
    medium <- medium_scenarios[as.character(medium_scenarios$medium_scenario_id) == scenario, , drop = FALSE]
    if (!nrow(medium) || (".no_constraints" %in% colnames(medium) && all(medium$.no_constraints))) medium <- NULL
    model <- rc_build_meta_module_gem(
      gem = gem,
      reaction_membership = reaction_membership,
      core_reactions = core_reactions,
      sample_id = "global",
      module_id = "GLOBAL_UNION",
      medium_table = medium,
      condition = NULL,
      target_direction = target_direction,
      solver = solver,
      time_limit = time_limit,
      fastcore_epsilon = fastcore_epsilon,
      max_support_reactions = max_support_reactions,
      strict = strict
    )
    model$shared_across_units <- TRUE
    file <- file.path(cache_dir, paste0("global_meta_module__medium_", safe(scenario), ".rds"))
    saveRDS(model, file)
    summaries[[scenario]] <- data.frame(
      medium_scenario = scenario,
      file = file,
      n_reactions = ncol(model$S),
      n_metabolites = nrow(model$S),
      n_biological_reactions = model$build_params$n_biological_reactions,
      n_fastcore_support_reactions = model$build_params$n_fastcore_support_reactions,
      target_status = model$target_status,
      build_strategy = "global_meta_module_gem",
      stringsAsFactors = FALSE
    )
    if (!nrow(model$target_directions)) next
    for (i in seq_len(nrow(model$target_directions))) {
      reaction <- as.character(model$target_directions$reaction_id[[i]])
      direction <- as.character(model$target_directions$target_direction[[i]])
      key <- paste0(
        "reaction=", utils::URLencode(reaction, reserved = TRUE),
        "::direction=", direction,
        "::medium=", utils::URLencode(scenario, reserved = TRUE)
      )
      cache[[key]] <- list(
        sample_id = "global",
        module_id = "GLOBAL_UNION",
        reaction_id = reaction,
        target_direction = direction,
        medium_scenario = scenario,
        condition = "all",
        file = file,
        build_strategy = "global_meta_module_gem"
      )
    }
  }
  attr(cache, "summary") <- .rc_bind_frames_fill(summaries)
  cache
}

#' Run COMPASS-like directional LPs with one shared structural GEM
#' @export
rc_run_microcompass <- function(layer1, gem, target_reactions = NULL,
                                 medium_table = NULL, medium_scenarios = NULL,
                                 mode = c("full_gem", "meta_module_gem"),
                                 reaction_membership = NULL, core_reactions = NULL,
                                 unit = c("metacell", "sample_celltype"),
                                 condition_col = "condition", sample_col = "sample_id",
                                 celltype_col = "cell_type", model_params = list(),
                                 penalty_weights = c(expr = 1.0, confidence = 0.5, missing = 1.0),
                                 omega = 0.95,
                                 target_direction = c("both", "forward", "reverse"),
                                 parallel = TRUE,
                                 solver = c("highs", "gurobi", "glpk"),
                                 time_limit = 60, flux_threshold = 1e-8,
                                 BPPARAM = NULL) {
  mode <- match.arg(mode)
  unit <- match.arg(unit)
  solver <- match.arg(solver)
  target_direction <- match.arg(target_direction)
  medium_scenarios <- .rc_validate_shared_medium(medium_scenarios %||% medium_table)
  matrices <- rc_layer2_unit_matrices(
    layer1,
    if (identical(unit, "metacell")) "metacell" else "sample_celltype",
    sample_col, celltype_col, condition_col
  )
  gem <- rc_annotate_reaction_roles(gem)
  direction_diagnostics <- NULL

  if (identical(mode, "full_gem")) {
    if (is.null(target_reactions) || !length(target_reactions)) stop("`target_reactions` is required in full-GEM mode.", call. = FALSE)
    directions <- rc_prepare_directional_targets(gem, target_reactions, target_direction)
    direction_diagnostics <- directions
    directions <- directions[directions$target_direction %in% c("forward", "reverse"), , drop = FALSE]
    if (!nrow(directions)) stop("No target reaction directions are allowed by the GEM bounds.", call. = FALSE)
    model_cache <- rc_build_full_gem_cache(
      gem = gem,
      dirs = directions,
      medium_scenarios = medium_scenarios,
      cache_dir = model_params$cache_dir %||% tempfile("RegCompassR_full_gem_cache_"),
      conditions = "all"
    )
  } else {
    model_cache <- .rc_build_global_meta_module_gem_cache(
      gem = gem,
      reaction_membership = reaction_membership,
      core_reactions = core_reactions,
      target_reactions = target_reactions,
      medium_scenarios = medium_scenarios,
      cache_dir = model_params$cache_dir %||% tempfile("RegCompassR_global_meta_module_cache_"),
      target_direction = target_direction,
      solver = solver,
      time_limit = model_params$completion_time_limit %||% max(time_limit, 300),
      fastcore_epsilon = model_params$fastcore_epsilon %||% 1e-4,
      max_support_reactions = model_params$max_support_reactions %||% 2000,
      strict = model_params$strict %||% TRUE
    )
    if (!length(model_cache)) stop("No parent-feasible global meta-module targets were available.", call. = FALSE)
    directions <- unique(do.call(rbind, lapply(model_cache, function(entry) {
      data.frame(reaction_id = entry$reaction_id, target_direction = entry$target_direction,
                 medium_scenario = entry$medium_scenario, stringsAsFactors = FALSE)
    })))
  }

  all_reactions <- unique(unlist(lapply(model_cache, function(entry) colnames(.rc_cache_gem(entry)$S)), use.names = FALSE))
  penalties <- rc_compute_multiome_penalty(
    rc_align_layer2_evidence(matrices$C_rel, all_reactions, NA_real_),
    rc_align_layer2_evidence(matrices$reaction_confidence, all_reactions, NA_real_),
    layer1$gpr_diagnostics,
    gem$reaction_roles,
    weights = penalty_weights
  )
  units <- colnames(matrices$C_rel)
  row_ids <- names(model_cache)
  penalty <- vmax <- matrix(NA_real_, nrow = length(row_ids), ncol = length(units), dimnames = list(row_ids, units))
  feasible <- evaluated <- matrix(FALSE, nrow = length(row_ids), ncol = length(units), dimnames = list(row_ids, units))
  tasks <- expand.grid(row_id = row_ids, unit_id = units, stringsAsFactors = FALSE)

  run_one <- function(task) {
    entry <- model_cache[[as.character(task$row_id)]]
    model <- .rc_cache_gem(entry)
    unit_id <- as.character(task$unit_id)
    answer <- rc_compass_two_step_lp_directional(
      S = model$S, lb = model$lb, ub = model$ub,
      target_reaction = entry$reaction_id,
      penalties = penalties$penalty[colnames(model$S), unit_id],
      target_direction = entry$target_direction,
      omega = omega, solver = solver, time_limit = time_limit,
      flux_threshold = flux_threshold
    )
    list(
      row_id = as.character(task$row_id),
      unit_id = unit_id,
      penalty = answer$penalty,
      vmax = answer$vmax,
      feasible = isTRUE(answer$feasible),
      diagnostics = data.frame(
        row_id = as.character(task$row_id), unit_id = unit_id,
        sample_id = "global", module_id = if (identical(mode, "meta_module_gem")) "GLOBAL_UNION" else NA_character_,
        reaction_id = entry$reaction_id, target_direction = entry$target_direction,
        medium_scenario = entry$medium_scenario, condition = "all",
        strict_feasible = isTRUE(answer$feasible), solver_status = answer$solver_status,
        step1_status = answer$step1_status, step2_status = answer$step2_status,
        target_status = model$target_status %||% if (isTRUE(answer$feasible)) "ok" else "structurally_infeasible",
        objective_value = answer$penalty, vmax = answer$vmax,
        stringsAsFactors = FALSE
      )
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
  lp_diagnostics <- do.call(rbind, lapply(results, `[[`, "diagnostics"))
  model_files <- unique(vapply(model_cache, function(entry) entry$file %||% "", character(1)))
  model_files <- model_files[nzchar(model_files)]
  model_diagnostics <- .rc_bind_frames_fill(lapply(model_files, function(file) {
    model <- readRDS(file)
    model$closure_diagnostics %||% data.frame()
  }))

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
      shared_gem = TRUE,
      shared_gem_scope = "all_metacells",
      flux_threshold = flux_threshold
    ),
    method = if (identical(mode, "full_gem")) "microCOMPASS shared full-GEM directional LP" else "microCOMPASS shared global meta-module-GEM directional LP"
  )
}

#' Summarize a microCOMPASS result
#' @export
rc_summarize_microcompass <- function(result) {
  evaluated <- result$evaluated %||% matrix(TRUE, nrow = nrow(result$feasible), ncol = ncol(result$feasible), dimnames = dimnames(result$feasible))
  data.frame(
    model_mode = result$model_mode %||% NA_character_,
    n_targets = nrow(result$score),
    n_units = ncol(result$score),
    n_evaluated = sum(evaluated),
    feasible_fraction = if (any(evaluated)) mean(result$feasible[evaluated]) else NA_real_,
    shared_gem = isTRUE(result$params$shared_gem),
    stringsAsFactors = FALSE
  )
}
