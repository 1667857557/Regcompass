#' Prepare directional targets for reversible reactions
#' @export
rc_prepare_directional_targets <- function(gem, target_reactions, mode = c("both_for_reversible", "forward_only", "user_defined")) {
  mode <- match.arg(mode); gv <- rc_validate_gem(gem); meta <- gem$reaction_meta
  rev <- stats::setNames(gv$lb < 0 & gv$ub > 0, gv$reactions)
  if (!is.null(meta) && "reversible" %in% colnames(meta)) rev[as.character(meta$reaction_id)] <- as.logical(meta$reversible)
  do.call(rbind, lapply(intersect(target_reactions, gv$reactions), function(r) {
    dirs <- if (mode == "both_for_reversible" && isTRUE(rev[[r]])) c("forward", "reverse") else "forward"
    data.frame(reaction_id = r, target_direction = dirs, reversible = isTRUE(rev[[r]]), stringsAsFactors = FALSE)
  }))
}
rc_compass_vmax_directional <- function(S, lb, ub, target_reaction, direction = "forward", solver = "highs", time_limit = 60) {
  j <- match(target_reaction, colnames(S)); c1 <- rep(0, ncol(S)); c1[j] <- if (direction == "reverse") 1 else -1
  step <- rc_solve_lp(c1, S, rep(0, nrow(S)), rep(0, nrow(S)), lb, ub, solver, time_limit)
  vmax <- if (identical(step$status, "optimal")) -step$objective_value else NA_real_
  list(feasible = is.finite(vmax) && vmax > 1e-8, vmax = vmax, status = step$status, runtime = step$runtime)
}
rc_compass_two_step_lp_directional <- function(S, lb, ub, target_reaction, penalties, target_direction = "forward", omega = 0.95, solver = "highs", time_limit = 60) {
  v <- rc_compass_vmax_directional(S, lb, ub, target_reaction, target_direction, solver, time_limit)
  if (!isTRUE(v$feasible)) return(list(feasible = FALSE, vmax = v$vmax, penalty = NA_real_, solver_status = v$status, step1_status = v$status, step2_status = "not_run", number_constraints = nrow(S), number_variables = ncol(S), runtime = v$runtime))
  j <- match(target_reaction, colnames(S)); step2 <- rc_build_abs_penalty_lp(S, lb, ub, penalties, j, omega * v$vmax)
  if (target_direction == "reverse") { A <- rep(0, 2*ncol(S)); A[j] <- -1; A[ncol(S)+j] <- 1; step2$A[nrow(step2$A),] <- A }
  ans <- rc_solve_lp(step2$obj, step2$A, step2$lhs, step2$rhs, step2$lb, step2$ub, solver, time_limit)
  list(feasible = identical(ans$status, "optimal"), vmax = v$vmax, penalty = if (identical(ans$status, "optimal")) ans$objective else NA_real_, solver_status = ans$status, step1_status = v$status, step2_status = ans$status, number_constraints = nrow(step2$A), number_variables = length(step2$obj), runtime = sum(c(v$runtime, ans$runtime), na.rm = TRUE))
}
#' Run target-local microCOMPASS analyses
#' @export
rc_run_microcompass <- function(layer1, gem, target_reactions,
                                medium_table = NULL, medium_scenarios = NULL,
                                unit = c("sample_celltype", "metacell"),
                                condition_col = "condition",
                                sample_col = "sample_id",
                                celltype_col = "cell_type",
                                microgem_params = list(),
                                penalty_weights = c(expr = 1.0, confidence = 0.5, missing = 1.0),
                                omega = 0.95,
                                target_direction = "both",
                                parallel = TRUE,
                                solver = c("highs", "gurobi", "glpk"),
                                time_limit = 60,
                                BPPARAM = NULL) {
  unit <- match.arg(unit)
  solver <- match.arg(solver)
  mats <- rc_layer2_unit_matrices(layer1, if (unit == "metacell") "metacell" else "sample_celltype", sample_col, celltype_col, condition_col)
  gem <- rc_annotate_reaction_roles(gem)
  dirs <- rc_prepare_directional_targets(gem, target_reactions, if (target_direction == "both") "both_for_reversible" else "forward_only")
  if (is.null(medium_scenarios)) {
    medium_scenarios <- medium_table
    if (is.null(medium_scenarios)) {
      medium_scenarios <- data.frame(medium_scenario_id = "base", exchange_reaction_id = character(), lb = numeric(), ub = numeric(), available = logical(), stringsAsFactors = FALSE)
    }
    if (!"medium_scenario_id" %in% colnames(medium_scenarios)) medium_scenarios$medium_scenario_id <- "custom"
  }
  cache_strategy <- microgem_params$strategy %||% "target_khop"
  if (identical(cache_strategy, "module_meso_gem")) {
    module_params <- microgem_params
    module_params$strategy <- NULL
    cache_dir <- module_params$cache_dir %||% tempfile("RegCompassR_module_gem_cache_")
    module_col <- module_params$module_col %||% "metabolic_module"
    force_cache <- module_params$force_cache %||% FALSE
    module_params$cache_dir <- NULL
    module_params$module_col <- NULL
    module_params$force_cache <- NULL
    mg_cache <- rc_build_module_gem_cache(
      gem = gem,
      dirs = dirs,
      medium_scenarios = medium_scenarios,
      cache_dir = cache_dir,
      module_col = module_col,
      module_gem_params = module_params,
      force = force_cache
    )
  } else {
    target_params <- microgem_params
    target_params$strategy <- NULL
    mg_cache <- rc_build_microgem_cache(gem = gem, dirs = dirs, medium_scenarios = medium_scenarios, microgem_params = target_params)
  }
  cache_gem <- function(entry) {
    if (is.list(entry) && !is.null(entry$file)) readRDS(entry$file) else entry
  }
  all_rxns <- unique(unlist(lapply(mg_cache, function(entry) colnames(cache_gem(entry)$S))))
  pen <- rc_compute_multiome_penalty(
    rc_align_layer2_evidence(mats$C_rel, all_rxns, NA_real_),
    rc_align_layer2_evidence(mats$reaction_confidence, all_rxns, NA_real_),
    layer1$gpr_diagnostics,
    gem$reaction_roles,
    weights = penalty_weights
  )
  units <- colnames(mats$C_rel)
  row_ids <- names(mg_cache)
  score <- penalty <- vmax <- matrix(NA_real_, nrow = length(row_ids), ncol = length(units), dimnames = list(row_ids, units))
  feasible <- matrix(FALSE, nrow = length(row_ids), ncol = length(units), dimnames = list(row_ids, units))
  tasks <- expand.grid(row_id = row_ids, unit_id = units, stringsAsFactors = FALSE)
  run_one <- function(task) {
    mg <- cache_gem(mg_cache[[task$row_id]])
    u <- task$unit_id
    parts <- strsplit(task$row_id, "::", fixed = TRUE)[[1]]
    reaction_id <- parts[[1]]
    target_dir <- parts[[2]]
    medium_scenario <- parts[[3]]
    p <- pen$penalty[colnames(mg$S), u]
    ans <- rc_compass_two_step_lp_directional(mg$S, mg$lb, mg$ub,
                                              target_reaction = reaction_id,
                                              penalties = p,
                                              target_direction = target_dir,
                                              omega = omega,
                                              solver = solver,
                                              time_limit = time_limit)
    list(row_id = task$row_id,
         unit_id = u,
         penalty = ans$penalty,
         vmax = ans$vmax,
         feasible = isTRUE(ans$feasible),
         diag = data.frame(reaction_id = reaction_id,
                           target_direction = target_dir,
                           medium_scenario = medium_scenario,
                           unit_id = u,
                           strict_feasible = isTRUE(ans$feasible),
                           solver_status = ans$solver_status,
                           step1_status = ans$step1_status,
                           step2_status = ans$step2_status,
                           objective_value = ans$penalty,
                           vmax = ans$vmax,
                           stringsAsFactors = FALSE))
  }
  res <- rc_parallel_lapply(split(tasks, seq_len(nrow(tasks))), function(x) run_one(x[1, ]), BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE)
  for (x in res) {
    penalty[x$row_id, x$unit_id] <- x$penalty
    vmax[x$row_id, x$unit_id] <- x$vmax
    feasible[x$row_id, x$unit_id] <- x$feasible
  }
  score <- rc_compass_score_from_penalty(penalty, feasible)
  lp_diagnostics <- do.call(rbind, lapply(res, `[[`, "diag"))
  microgem_diagnostics <- do.call(rbind, lapply(mg_cache, function(entry) cache_gem(entry)$closure_diagnostics))
  scenarios <- unique(as.character(medium_scenarios$medium_scenario_id))
  sens <- if (length(scenarios) > 1L) {
    aggregate(strict_feasible ~ reaction_id + target_direction + unit_id, lp_diagnostics, function(x) length(unique(x)) > 1L)
  } else {
    NULL
  }
  if (!is.null(sens)) names(sens)[names(sens) == "strict_feasible"] <- "medium_sensitive_flag"
  list(score = score,
       penalty = penalty,
       vmax = vmax,
       feasible = feasible,
       target_direction = dirs,
       medium_scenarios = medium_scenarios,
       medium_sensitivity_summary = sens,
       microgem_cache_summary = attr(mg_cache, "summary"),
       microgem_diagnostics = microgem_diagnostics,
       lp_diagnostics = lp_diagnostics,
       penalty_components = pen$components,
       evidence_policy = pen$evidence_policy,
       unit_meta = mats$unit_meta,
       params = list(unit = unit,
                     omega = omega,
                     target_direction = target_direction,
                     evidence_policy = "RNA+ATAC-GPR evidence affects penalties only, not structural micro-GEM membership.",
                     interpretation = "multiome-supported reaction capacity potential from strict cached microCOMPASS"),
       method = "microCOMPASS strict target-local reaction-potential LP")
}
#' Summarize microCOMPASS result
#' @export
rc_summarize_microcompass <- function(result) data.frame(n_targets = nrow(result$score), n_units = ncol(result$score), feasible_fraction = mean(result$feasible, na.rm = TRUE))
