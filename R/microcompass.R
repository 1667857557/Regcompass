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
  vmax <- if (identical(step$status, "optimal")) if (direction == "reverse") step$objective else -step$objective else NA_real_
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
rc_run_microcompass <- function(layer1, gem, target_reactions, medium_table = NULL, unit = c("sample_celltype", "metacell"), condition_col = "condition", sample_col = "sample_id", celltype_col = "cell_type", microgem_params = list(), penalty_weights = c(expr = 1.0, confidence = 0.5, missing = 1.0), omega = 0.95, target_direction = "both", run_relaxed = TRUE, run_fva = TRUE, solver = c("highs", "gurobi", "glpk"), time_limit = 60, BPPARAM = NULL) {
  unit <- match.arg(unit); solver <- match.arg(solver)
  mats <- rc_layer2_unit_matrices(layer1, if (unit == "metacell") "pool" else "sample_celltype", sample_col, celltype_col, condition_col)
  gem <- rc_annotate_reaction_roles(gem)
  dirs <- rc_prepare_directional_targets(gem, target_reactions, if (target_direction == "both") "both_for_reversible" else "forward_only")
  ids <- paste(dirs$reaction_id, dirs$target_direction, sep = "::"); units <- colnames(mats$C_rel)
  score <- penalty <- vmax <- matrix(NA_real_, length(ids), length(units), dimnames = list(ids, units)); feasible <- matrix(FALSE, length(ids), length(units), dimnames = list(ids, units))
  lpd <- list(); mgd <- list(); relaxed <- list(); fva <- list(); idx <- 0L
  all_rxns <- unique(unlist(lapply(dirs$reaction_id, function(r) colnames(do.call(rc_build_target_microgem, c(list(gem = gem, target_reaction = r, medium_table = medium_table), microgem_params))$S))))
  pen <- rc_compute_multiome_penalty(rc_align_layer2_evidence(mats$C_rel, all_rxns, NA_real_), rc_align_layer2_evidence(mats$reaction_confidence, all_rxns, NA_real_), layer1$gpr_diagnostics, gem$reaction_roles, weights = penalty_weights)
  for (i in seq_len(nrow(dirs))) {
    mg <- do.call(rc_build_target_microgem, c(list(gem = gem, target_reaction = dirs$reaction_id[i], medium_table = medium_table), microgem_params)); mgd[[i]] <- mg$closure_diagnostics
    for (u in units) {
      p <- pen$penalty[colnames(mg$S), u]; ans <- rc_compass_two_step_lp_directional(mg$S, mg$lb, mg$ub, dirs$reaction_id[i], p, dirs$target_direction[i], omega, solver, time_limit)
      row <- paste(dirs$reaction_id[i], dirs$target_direction[i], sep = "::"); penalty[row,u] <- ans$penalty; vmax[row,u] <- ans$vmax; feasible[row,u] <- isTRUE(ans$feasible); idx <- idx + 1L
      lpd[[idx]] <- data.frame(reaction_id=dirs$reaction_id[i], target_direction=dirs$target_direction[i], unit_id=u, solver_status=ans$solver_status, step1_status=ans$step1_status, step2_status=ans$step2_status, objective_value=ans$penalty, vmax=ans$vmax, stringsAsFactors=FALSE)
      if (run_relaxed) relaxed[[idx]] <- rc_run_relaxed_balance_lp(mg, p, dirs$reaction_id[i], dirs$target_direction[i], omega, solver = solver, time_limit = time_limit)
      if (run_fva) fva[[idx]] <- rc_run_selected_fva(mg, dirs$reaction_id[i], target_direction = dirs$target_direction[i], solver = solver, time_limit = time_limit)
    }
  }
  score[] <- rc_sigmoid(-(penalty - stats::median(penalty, na.rm = TRUE))); score[!feasible] <- 0
  list(score=score, penalty=penalty, vmax=vmax, feasible=feasible, target_direction=dirs, microgem_diagnostics=do.call(rbind, mgd), lp_diagnostics=do.call(rbind, lpd), penalty_components=pen$components, relaxed=if(run_relaxed) do.call(rbind, relaxed) else NULL, fva=if(run_fva) do.call(rbind, fva) else NULL, unit_meta=mats$unit_meta, params=list(unit=unit, omega=omega, target_direction=target_direction), method="microCOMPASS target-local two-step penalty LP")
}
#' @export
rc_run_layer2_compass_lp <- function(...) { .Deprecated("rc_run_microcompass"); rc_run_microcompass(...) }
#' Summarize microCOMPASS result
#' @export
rc_summarize_microcompass <- function(result) data.frame(n_targets = nrow(result$score), n_units = ncol(result$score), feasible_fraction = mean(result$feasible, na.rm = TRUE))
