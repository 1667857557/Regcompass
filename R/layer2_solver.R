rc_compass_two_step_lp <- function(S, lb, ub, target_reaction, penalties, omega = 0.95, solver = "highs", time_limit = 60, tol_vmax = 1e-8) {
  rxns <- colnames(S); j <- match(target_reaction, rxns)
  if (is.na(j)) stop("Target reaction is not in `S`.", call. = FALSE)
  c1 <- rep(0, length(rxns)); c1[j] <- -1
  step1 <- rc_solve_lp(c1, S, rep(0, nrow(S)), rep(0, nrow(S)), lb, ub, solver, time_limit)
  vmax <- if (identical(step1$status, "optimal")) -step1$objective else NA_real_
  if (!is.finite(vmax) || vmax <= tol_vmax) {
    return(list(feasible = FALSE, vmax = vmax, penalty = NA_real_, solver_status = step1$status,
                step1_status = step1$status, step2_status = "not_run",
                number_constraints = nrow(S), number_variables = ncol(S), runtime = step1$runtime %||% NA_real_))
  }
  step2 <- rc_build_abs_penalty_lp(S, lb, ub, penalties, target_index = j, target_min = omega * vmax)
  ans <- rc_solve_lp(step2$obj, step2$A, step2$lhs, step2$rhs, step2$lb, step2$ub, solver, time_limit)
  list(feasible = identical(ans$status, "optimal"), vmax = vmax, penalty = if (identical(ans$status, "optimal")) ans$objective else NA_real_,
       solver_status = ans$status, step1_status = step1$status, step2_status = ans$status,
       number_constraints = nrow(step2$A), number_variables = length(step2$obj),
       runtime = sum(c(step1$runtime, ans$runtime), na.rm = TRUE))
}

rc_build_abs_penalty_lp <- function(S, lb, ub, penalties, target_index, target_min) {
  n <- ncol(S); p <- as.numeric(penalties[colnames(S)]); p[!is.finite(p)] <- if (any(is.finite(p))) max(p[is.finite(p)], 20) else 20; p <- pmax(p, 0)
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

rc_standardize_lp_result <- function(status, objective_value, primal_solution, runtime, message = NULL) {
  map <- tolower(as.character(status %||% "error"))
  std <- if (grepl("optimal", map)) "optimal" else if (grepl("infeas", map)) "infeasible" else if (grepl("unbound", map)) "unbounded" else if (grepl("time|limit", map)) "time_limit" else if (identical(map, "optimal")) "optimal" else map
  structure(list(status = std,
                 objective_value = as.numeric(objective_value),
                 objective = as.numeric(objective_value),
                 primal_solution = primal_solution,
                 solution = primal_solution,
                 runtime = runtime,
                 message = message),
            objective_sense = "min")
}

.rc_highs_status <- function(res, n_variables) {
  text_candidates <- c(res$model_status, res$status, res$status_message, res$message)
  text_candidates <- as.character(unlist(text_candidates, use.names = FALSE))
  text_candidates <- text_candidates[!is.na(text_candidates) & nzchar(text_candidates)]
  text <- tolower(paste(text_candidates, collapse = " "))
  if (grepl("optimal", text)) return("optimal")
  if (grepl("infeasible", text)) return("infeasible")
  if (grepl("unbounded", text)) return("unbounded")
  if (grepl("time|limit", text)) return("time_limit")
  code <- suppressWarnings(as.integer(res$status)[1L])
  solution <- res$primal_solution %||% res$solution
  objective <- res$objective_value %||% res$objective
  if (identical(code, 7L) && length(objective) == 1L && is.finite(objective) && length(solution) == n_variables && all(is.finite(solution))) return("optimal")
  if (is.finite(code)) return(paste0("highs_status_", code))
  "error"
}

rc_solve_lp_gurobi <- function(obj, A, lhs, rhs, lb, ub, time_limit) {
  start_time <- proc.time()[["elapsed"]]
  if (!requireNamespace("gurobi", quietly = TRUE)) return(rc_standardize_lp_result("error", NA_real_, NULL, 0, "gurobi package not installed"))
  sense <- ifelse(is.finite(lhs) & is.finite(rhs) & lhs == rhs, "=", ifelse(is.finite(lhs), ">", "<"))
  b <- ifelse(sense == ">", lhs, rhs)
  model <- list(modelsense = "min", obj = obj, A = A, sense = sense, rhs = b, lb = lb, ub = ub)
  res <- tryCatch(gurobi::gurobi(model, params = list(TimeLimit = time_limit, OutputFlag = 0)), error = function(e) e)
  if (inherits(res, "error")) return(rc_standardize_lp_result("error", NA_real_, NULL, proc.time()[["elapsed"]] - start_time, conditionMessage(res)))
  rc_standardize_lp_result(res$status, res$objval, res$x, proc.time()[["elapsed"]] - start_time)
}

rc_solve_lp_highs <- function(obj, A, lhs, rhs, lb, ub, time_limit) {
  start_time <- proc.time()[["elapsed"]]
  if (!requireNamespace("highs", quietly = TRUE)) {
    return(rc_standardize_lp_result("error", NA_real_, NULL, 0, "highs package not installed"))
  }

  fmls <- names(formals(highs::highs_solve))
  args <- list(
    L = as.numeric(obj),
    lower = as.numeric(lb),
    upper = as.numeric(ub),
    A = A,
    maximum = FALSE
  )

  if ("Q" %in% fmls) args$Q <- NULL

  if (all(c("lhs", "rhs") %in% fmls)) {
    args$lhs <- as.numeric(lhs)
    args$rhs <- as.numeric(rhs)
  } else if (all(c("lower_bounds", "upper_bounds") %in% fmls)) {
    args$lower_bounds <- as.numeric(lb)
    args$upper_bounds <- as.numeric(ub)
    args$lower <- as.numeric(lhs)
    args$upper <- as.numeric(rhs)
  } else {
    return(rc_standardize_lp_result(
      "error",
      NA_real_,
      NULL,
      proc.time()[["elapsed"]] - start_time,
      "Unsupported highs::highs_solve() API; expected either lhs/rhs or lower_bounds/upper_bounds interface."
    ))
  }

  if ("types" %in% fmls) args$types <- rep("C", length(obj))
  if ("options" %in% fmls || "..." %in% fmls) args$options <- list(time_limit = time_limit)

  res <- tryCatch(do.call(highs::highs_solve, args), error = function(e) e)
  if (inherits(res, "error")) {
    return(rc_standardize_lp_result("error", NA_real_, NULL, proc.time()[["elapsed"]] - start_time, conditionMessage(res)))
  }
  status <- .rc_highs_status(res, n_variables = length(obj))
  out <- rc_standardize_lp_result(
    status,
    res$objective_value %||% res$objective,
    res$primal_solution %||% res$solution,
    proc.time()[["elapsed"]] - start_time,
    res$message %||% res$status_message %||% NULL
  )
  out$raw_status <- res$status %||% NA
  out$raw_model_status <- res$model_status %||% NA
  out$solver_version <- as.character(utils::packageVersion("highs"))
  out
}

rc_solve_lp_glpk <- function(obj, A, lhs, rhs, lb, ub, time_limit) {
  start_time <- proc.time()[["elapsed"]]
  if (!requireNamespace("Rglpk", quietly = TRUE)) return(rc_standardize_lp_result("error", NA_real_, NULL, 0, "Rglpk package not installed"))
  dir <- rep("==", length(lhs)); dir[is.infinite(lhs)] <- "<="; dir[is.infinite(rhs)] <- ">="; b <- ifelse(dir == ">=", lhs, rhs)
  res <- tryCatch(Rglpk::Rglpk_solve_LP(obj, A, dir, b, bounds = list(lower = list(ind = seq_along(lb), val = lb), upper = list(ind = seq_along(ub), val = ub)), max = FALSE, control = list(tm_limit = time_limit * 1000)), error = function(e) e)
  if (inherits(res, "error")) return(rc_standardize_lp_result("error", NA_real_, NULL, proc.time()[["elapsed"]] - start_time, conditionMessage(res)))
  status <- if (res$status == 0) "optimal" else paste0("glpk_status_", res$status)
  rc_standardize_lp_result(status, res$optimum, res$solution, proc.time()[["elapsed"]] - start_time)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
