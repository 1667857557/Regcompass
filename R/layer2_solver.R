# Directional COMPASS and generic LP solver helpers.

rc_compass_vmax_directional <- function(S, lb, ub, target_reaction,
                                         direction = c("forward", "reverse"),
                                         solver = "highs", time_limit = 60,
                                         flux_threshold = 1e-8) {
  direction <- match.arg(direction)
  if (!is.finite(flux_threshold) || flux_threshold <= 0) {
    stop("`flux_threshold` must be positive.", call. = FALSE)
  }
  rxns <- colnames(S)
  j <- match(target_reaction, rxns)
  if (is.na(j)) stop("Target reaction is not in `S`.", call. = FALSE)

  lb <- rc_align_bound(lb, rxns, default = -1000, name = "lb")
  ub <- rc_align_bound(ub, rxns, default = 1000, name = "ub")
  direction_sign <- if (identical(direction, "reverse")) -1 else 1
  direction_allowed <- if (direction_sign > 0) {
    ub[[j]] >= flux_threshold
  } else {
    lb[[j]] <= -flux_threshold
  }
  if (!isTRUE(direction_allowed)) {
    return(list(
      feasible = FALSE, vmax = 0,
      status = "direction_not_allowed", runtime = 0
    ))
  }

  objective <- rep(0, length(rxns))
  objective[[j]] <- -direction_sign
  step <- rc_solve_lp(
    obj = objective,
    A = S,
    lhs = rep(0, nrow(S)),
    rhs = rep(0, nrow(S)),
    lb = lb,
    ub = ub,
    solver = solver,
    time_limit = time_limit
  )
  vmax <- if (identical(step$status, "optimal")) {
    -step$objective_value
  } else {
    NA_real_
  }
  list(
    feasible = is.finite(vmax) && vmax >= flux_threshold * (1 - 1e-7),
    vmax = vmax,
    status = step$status,
    runtime = step$runtime %||% NA_real_
  )
}

rc_build_abs_penalty_lp <- function(S, lb, ub, penalties, target_index,
                                    target_min,
                                    target_direction = c("forward", "reverse"),
                                    penalty_floor = 1e-12) {
  target_direction <- match.arg(target_direction)
  n <- ncol(S)
  rxns <- colnames(S)
  if (!length(target_index) || is.na(target_index) ||
      target_index < 1L || target_index > n) {
    stop("`target_index` is invalid.", call. = FALSE)
  }
  lb <- rc_align_bound(lb, rxns, default = -1000, name = "lb")
  ub <- rc_align_bound(ub, rxns, default = 1000, name = "ub")

  if (is.null(names(penalties)) && length(penalties) == n) {
    p <- as.numeric(penalties)
  } else {
    p <- as.numeric(penalties[rxns])
  }
  finite_p <- p[is.finite(p)]
  replacement <- if (length(finite_p)) max(c(finite_p, 20)) else 20
  p[!is.finite(p)] <- replacement
  p <- pmax(p, penalty_floor)

  S <- methods::as(S, "dgCMatrix")
  identity <- Matrix::Diagonal(n)
  zero <- Matrix::Matrix(0, nrow = nrow(S), ncol = n, sparse = TRUE)
  mass_balance <- cbind(S, zero)
  abs_positive <- cbind(identity, -identity)
  abs_negative <- cbind(-identity, -identity)
  target <- Matrix::Matrix(0, nrow = 1L, ncol = 2L * n, sparse = TRUE)
  target[1L, target_index] <- if (identical(target_direction, "reverse")) -1 else 1

  A <- rbind(mass_balance, abs_positive, abs_negative, target)
  lhs <- c(rep(0, nrow(S)), rep(-Inf, 2L * n), target_min)
  rhs <- c(rep(0, nrow(S)), rep(0, 2L * n), Inf)
  absolute_upper <- pmax(abs(lb), abs(ub))

  list(
    obj = c(rep(0, n), p),
    A = A,
    lhs = lhs,
    rhs = rhs,
    lb = c(lb, rep(0, n)),
    ub = c(ub, absolute_upper),
    n_flux_variables = n
  )
}

rc_compass_two_step_lp_directional <- function(S, lb, ub, target_reaction,
                                                penalties,
                                                target_direction = c("forward", "reverse"),
                                                omega = 0.95,
                                                solver = "highs",
                                                time_limit = 60,
                                                flux_threshold = 1e-8) {
  target_direction <- match.arg(target_direction)
  if (!is.finite(omega) || omega <= 0 || omega > 1) {
    stop("`omega` must lie in (0, 1].", call. = FALSE)
  }
  vmax_result <- rc_compass_vmax_directional(
    S = S,
    lb = lb,
    ub = ub,
    target_reaction = target_reaction,
    direction = target_direction,
    solver = solver,
    time_limit = time_limit,
    flux_threshold = flux_threshold
  )
  if (!isTRUE(vmax_result$feasible)) {
    return(list(
      feasible = FALSE,
      vmax = vmax_result$vmax,
      penalty = NA_real_,
      solver_status = vmax_result$status,
      step1_status = vmax_result$status,
      step2_status = "not_run",
      number_constraints = nrow(S),
      number_variables = ncol(S),
      runtime = vmax_result$runtime
    ))
  }

  target_index <- match(target_reaction, colnames(S))
  lp <- rc_build_abs_penalty_lp(
    S = S,
    lb = lb,
    ub = ub,
    penalties = penalties,
    target_index = target_index,
    target_min = omega * vmax_result$vmax,
    target_direction = target_direction
  )
  answer <- rc_solve_lp(
    obj = lp$obj,
    A = lp$A,
    lhs = lp$lhs,
    rhs = lp$rhs,
    lb = lp$lb,
    ub = lp$ub,
    solver = solver,
    time_limit = time_limit
  )

  list(
    feasible = identical(answer$status, "optimal"),
    vmax = vmax_result$vmax,
    penalty = if (identical(answer$status, "optimal")) answer$objective else NA_real_,
    solver_status = answer$status,
    step1_status = vmax_result$status,
    step2_status = answer$status,
    number_constraints = nrow(lp$A),
    number_variables = length(lp$obj),
    runtime = sum(c(vmax_result$runtime, answer$runtime), na.rm = TRUE)
  )
}

rc_compass_two_step_lp <- function(S, lb, ub, target_reaction, penalties,
                                   omega = 0.95, solver = "highs",
                                   time_limit = 60, tol_vmax = 1e-8) {
  rc_compass_two_step_lp_directional(
    S = S,
    lb = lb,
    ub = ub,
    target_reaction = target_reaction,
    penalties = penalties,
    target_direction = "forward",
    omega = omega,
    solver = solver,
    time_limit = time_limit,
    flux_threshold = tol_vmax
  )
}

rc_solve_lp <- function(obj, A, lhs, rhs, lb, ub,
                        solver = "highs", time_limit = 60) {
  if (solver == "gurobi") {
    return(rc_solve_lp_gurobi(obj, A, lhs, rhs, lb, ub, time_limit))
  }
  if (solver == "highs") {
    return(rc_solve_lp_highs(obj, A, lhs, rhs, lb, ub, time_limit))
  }
  if (solver == "glpk") {
    return(rc_solve_lp_glpk(obj, A, lhs, rhs, lb, ub, time_limit))
  }
  stop("Unsupported solver.", call. = FALSE)
}

rc_standardize_lp_result <- function(status, objective_value,
                                     primal_solution, runtime,
                                     message = NULL) {
  status_text <- tolower(as.character(status %||% "error"))
  standardized <- if (grepl("optimal", status_text)) {
    "optimal"
  } else if (grepl("infeas", status_text)) {
    "infeasible"
  } else if (grepl("unbound", status_text)) {
    "unbounded"
  } else if (grepl("time|limit", status_text)) {
    "time_limit"
  } else {
    status_text
  }
  structure(
    list(
      status = standardized,
      objective_value = as.numeric(objective_value),
      objective = as.numeric(objective_value),
      primal_solution = primal_solution,
      solution = primal_solution,
      runtime = runtime,
      message = message
    ),
    objective_sense = "min"
  )
}

.rc_highs_status <- function(result, n_variables) {
  text_candidates <- c(
    result$model_status, result$status,
    result$status_message, result$message
  )
  text_candidates <- as.character(unlist(text_candidates, use.names = FALSE))
  text_candidates <- text_candidates[
    !is.na(text_candidates) & nzchar(text_candidates)
  ]
  text <- tolower(paste(text_candidates, collapse = " "))
  if (grepl("optimal", text)) return("optimal")
  if (grepl("infeasible", text)) return("infeasible")
  if (grepl("unbounded", text)) return("unbounded")
  if (grepl("time|limit", text)) return("time_limit")

  code <- suppressWarnings(as.integer(result$status)[1L])
  solution <- result$primal_solution %||% result$solution
  objective <- result$objective_value %||% result$objective
  if (identical(code, 7L) && length(objective) == 1L &&
      is.finite(objective) && length(solution) == n_variables &&
      all(is.finite(solution))) {
    return("optimal")
  }
  if (is.finite(code)) return(paste0("highs_status_", code))
  "error"
}

rc_solve_lp_gurobi <- function(obj, A, lhs, rhs, lb, ub, time_limit) {
  started <- proc.time()[["elapsed"]]
  if (!requireNamespace("gurobi", quietly = TRUE)) {
    return(rc_standardize_lp_result(
      "error", NA_real_, NULL, 0,
      "gurobi package not installed"
    ))
  }
  sense <- ifelse(
    is.finite(lhs) & is.finite(rhs) & lhs == rhs,
    "=",
    ifelse(is.finite(lhs), ">", "<")
  )
  bound <- ifelse(sense == ">", lhs, rhs)
  model <- list(
    modelsense = "min", obj = obj, A = A,
    sense = sense, rhs = bound, lb = lb, ub = ub
  )
  result <- tryCatch(
    gurobi::gurobi(
      model,
      params = list(TimeLimit = time_limit, OutputFlag = 0)
    ),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    return(rc_standardize_lp_result(
      "error", NA_real_, NULL,
      proc.time()[["elapsed"]] - started,
      conditionMessage(result)
    ))
  }
  rc_standardize_lp_result(
    result$status, result$objval, result$x,
    proc.time()[["elapsed"]] - started
  )
}

rc_solve_lp_highs <- function(obj, A, lhs, rhs, lb, ub, time_limit) {
  started <- proc.time()[["elapsed"]]
  if (!requireNamespace("highs", quietly = TRUE)) {
    return(rc_standardize_lp_result(
      "error", NA_real_, NULL, 0,
      "highs package not installed"
    ))
  }

  formals_highs <- names(formals(highs::highs_solve))
  arguments <- list(
    L = as.numeric(obj),
    lower = as.numeric(lb),
    upper = as.numeric(ub),
    A = A,
    maximum = FALSE
  )
  if ("Q" %in% formals_highs) arguments$Q <- NULL
  if (all(c("lhs", "rhs") %in% formals_highs)) {
    arguments$lhs <- as.numeric(lhs)
    arguments$rhs <- as.numeric(rhs)
  } else if (all(c("lower_bounds", "upper_bounds") %in% formals_highs)) {
    arguments$lower_bounds <- as.numeric(lb)
    arguments$upper_bounds <- as.numeric(ub)
    arguments$lower <- as.numeric(lhs)
    arguments$upper <- as.numeric(rhs)
  } else {
    return(rc_standardize_lp_result(
      "error", NA_real_, NULL,
      proc.time()[["elapsed"]] - started,
      paste(
        "Unsupported highs::highs_solve() API; expected",
        "either lhs/rhs or lower_bounds/upper_bounds interface."
      )
    ))
  }
  if ("types" %in% formals_highs) {
    arguments$types <- rep("C", length(obj))
  }
  if ("options" %in% formals_highs || "..." %in% formals_highs) {
    arguments$options <- list(time_limit = time_limit)
  }

  result <- tryCatch(
    do.call(highs::highs_solve, arguments),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    return(rc_standardize_lp_result(
      "error", NA_real_, NULL,
      proc.time()[["elapsed"]] - started,
      conditionMessage(result)
    ))
  }
  status <- .rc_highs_status(result, n_variables = length(obj))
  output <- rc_standardize_lp_result(
    status,
    result$objective_value %||% result$objective,
    result$primal_solution %||% result$solution,
    proc.time()[["elapsed"]] - started,
    result$message %||% result$status_message %||% NULL
  )
  output$raw_status <- result$status %||% NA
  output$raw_model_status <- result$model_status %||% NA
  output$solver_version <- as.character(utils::packageVersion("highs"))
  output
}

rc_solve_lp_glpk <- function(obj, A, lhs, rhs, lb, ub, time_limit) {
  started <- proc.time()[["elapsed"]]
  if (!requireNamespace("Rglpk", quietly = TRUE)) {
    return(rc_standardize_lp_result(
      "error", NA_real_, NULL, 0,
      "Rglpk package not installed"
    ))
  }
  equality <- is.finite(lhs) & is.finite(rhs) & abs(lhs - rhs) <= 1e-12
  lower_only <- is.finite(lhs) & is.infinite(rhs)
  upper_only <- is.infinite(lhs) & is.finite(rhs)
  if (any(!(equality | lower_only | upper_only))) {
    return(rc_standardize_lp_result(
      "error", NA_real_, NULL, 0,
      "Rglpk adapter only supports equality and one-sided rows"
    ))
  }
  direction <- ifelse(equality, "==", ifelse(lower_only, ">=", "<="))
  bound <- ifelse(lower_only, lhs, rhs)
  result <- tryCatch(
    Rglpk::Rglpk_solve_LP(
      obj, A, direction, bound,
      bounds = list(
        lower = list(ind = seq_along(lb), val = lb),
        upper = list(ind = seq_along(ub), val = ub)
      ),
      max = FALSE,
      control = list(tm_limit = time_limit * 1000)
    ),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    return(rc_standardize_lp_result(
      "error", NA_real_, NULL,
      proc.time()[["elapsed"]] - started,
      conditionMessage(result)
    ))
  }
  status <- if (result$status == 0) {
    "optimal"
  } else {
    paste0("glpk_status_", result$status)
  }
  rc_standardize_lp_result(
    status, result$optimum, result$solution,
    proc.time()[["elapsed"]] - started
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x
