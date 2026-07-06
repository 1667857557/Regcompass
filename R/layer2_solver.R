rc_compass_two_step_lp <- function(S, lb, ub, target_reaction, penalties, omega = 0.95, solver = "highs", time_limit = 60, tol_vmax = 1e-8) {
  rxns <- colnames(S); j <- match(target_reaction, rxns)
  if (is.na(j)) stop("Target reaction is not in `S`.", call. = FALSE)
  c1 <- rep(0, length(rxns)); c1[j] <- -1
  step1 <- rc_solve_lp(c1, S, rep(0, nrow(S)), rep(0, nrow(S)), lb, ub, solver, time_limit)
  vmax <- if (identical(step1$status, "optimal")) -step1$objective else NA_real_
  if (!is.finite(vmax) || vmax <= tol_vmax) return(list(feasible = FALSE, vmax = vmax, penalty = NA_real_, solver_status = step1$status))
  step2 <- rc_build_abs_penalty_lp(S, lb, ub, penalties, target_index = j, target_min = omega * vmax)
  ans <- rc_solve_lp(step2$obj, step2$A, step2$lhs, step2$rhs, step2$lb, step2$ub, solver, time_limit)
  list(feasible = identical(ans$status, "optimal"), vmax = vmax, penalty = if (identical(ans$status, "optimal")) ans$objective else NA_real_, solver_status = ans$status)
}

rc_build_abs_penalty_lp <- function(S, lb, ub, penalties, target_index, target_min) {
  n <- ncol(S); p <- as.numeric(penalties[colnames(S)]); p[!is.finite(p)] <- max(p[is.finite(p)], 20)
  # variables are vplus, vminus; v = vplus - vminus
  Aeq <- cbind(S, -S); lhs <- rhs <- rep(0, nrow(S))
  A_target <- rep(0, 2*n); A_target[target_index] <- 1; A_target[n + target_index] <- -1
  A <- rbind(Aeq, A_target); lhs <- c(lhs, target_min); rhs <- c(rhs, Inf)
  list(obj = c(p, p), A = A, lhs = lhs, rhs = rhs, lb = rep(0, 2 * n), ub = c(pmax(ub, 0), pmax(-lb, 0)))
}

rc_solve_lp <- function(obj, A, lhs, rhs, lb, ub, solver = "highs", time_limit = 60) {
  if (solver == "gurobi") return(rc_solve_lp_gurobi(obj, A, lhs, rhs, lb, ub, time_limit))
  if (solver == "highs") return(rc_solve_lp_highs(obj, A, lhs, rhs, lb, ub, time_limit))
  if (solver == "glpk") return(rc_solve_lp_glpk(obj, A, lhs, rhs, lb, ub, time_limit))
  stop("Unsupported solver.", call. = FALSE)
}

rc_solve_lp_gurobi <- function(obj, A, lhs, rhs, lb, ub, time_limit) {
  if (!requireNamespace("gurobi", quietly = TRUE)) return(list(status = "error", objective = NA_real_, solution = NULL, message = "gurobi package not installed"))
  sense <- ifelse(is.finite(lhs) & is.finite(rhs) & lhs == rhs, "=", ifelse(is.finite(lhs), ">", "<"))
  b <- ifelse(sense == ">", lhs, rhs)
  model <- list(modelsense = "min", obj = obj, A = A, sense = sense, rhs = b, lb = lb, ub = ub)
  res <- tryCatch(gurobi::gurobi(model, params = list(TimeLimit = time_limit, OutputFlag = 0)), error = function(e) e)
  if (inherits(res, "error")) return(list(status = "error", objective = NA_real_, solution = NULL, message = conditionMessage(res)))
  status <- if (tolower(res$status) == "optimal") "optimal" else tolower(res$status)
  list(status = status, objective = res$objval, solution = res$x)
}

rc_solve_lp_highs <- function(obj, A, lhs, rhs, lb, ub, time_limit) {
  if (!requireNamespace("highs", quietly = TRUE)) return(list(status = "error", objective = NA_real_, solution = NULL, message = "highs package not installed"))
  res <- tryCatch(highs::highs_solve(L = obj, lower = lhs, upper = rhs, A = A, lower_bounds = lb, upper_bounds = ub, maximum = FALSE), error = function(e) e)
  if (inherits(res, "error")) return(list(status = "error", objective = NA_real_, solution = NULL, message = conditionMessage(res)))
  status <- if (!is.null(res$status) && grepl("optimal", tolower(res$status))) "optimal" else tolower(as.character(res$status %||% "error"))
  list(status = status, objective = as.numeric(res$objective_value %||% res$objective), solution = res$primal_solution %||% res$solution)
}

rc_solve_lp_glpk <- function(obj, A, lhs, rhs, lb, ub, time_limit) {
  if (!requireNamespace("Rglpk", quietly = TRUE)) return(list(status = "error", objective = NA_real_, solution = NULL, message = "Rglpk package not installed"))
  dir <- rep("==", length(lhs)); dir[is.infinite(lhs)] <- "<="; dir[is.infinite(rhs)] <- ">="; b <- ifelse(dir == ">=", lhs, rhs)
  res <- tryCatch(Rglpk::Rglpk_solve_LP(obj, A, dir, b, bounds = list(lower = list(ind = seq_along(lb), val = lb), upper = list(ind = seq_along(ub), val = ub)), max = FALSE, control = list(tm_limit = time_limit * 1000)), error = function(e) e)
  if (inherits(res, "error")) return(list(status = "error", objective = NA_real_, solution = NULL, message = conditionMessage(res)))
  status <- if (res$status == 0) "optimal" else paste0("glpk_status_", res$status)
  list(status = status, objective = res$optimum, solution = res$solution)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
