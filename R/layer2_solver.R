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
  S <- .rc_as_dgCMatrix(S)
  n <- ncol(S)
  reactions <- colnames(S)
  if (is.null(reactions) || anyNA(reactions) || any(!nzchar(reactions)) ||
      anyDuplicated(reactions)) {
    stop("`S` must have unique non-empty reaction IDs in colnames().", call. = FALSE)
  }
  if (length(target_index) != 1L || is.na(target_index) ||
      target_index < 1L || target_index > n) {
    stop("`target_index` is invalid.", call. = FALSE)
  }
  if (!is.numeric(target_min) || length(target_min) != 1L ||
      !is.finite(target_min) || target_min < 0) {
    stop("`target_min` must be one finite non-negative number.", call. = FALSE)
  }
  if (!is.numeric(penalty_floor) || length(penalty_floor) != 1L ||
      !is.finite(penalty_floor) || penalty_floor < 0) {
    stop("`penalty_floor` must be one finite non-negative number.", call. = FALSE)
  }
  lb <- rc_align_bound(lb, reactions, default = -1000, name = "lb")
  ub <- rc_align_bound(ub, reactions, default = 1000, name = "ub")

  if (is.null(names(penalties))) {
    if (length(penalties) != n) {
      stop("Unnamed `penalties` must have one value per reaction.", call. = FALSE)
    }
    penalty <- as.numeric(penalties)
  } else {
    penalty_names <- as.character(names(penalties))
    if (anyNA(penalty_names) || any(!nzchar(penalty_names)) ||
        anyDuplicated(penalty_names)) {
      stop("Named `penalties` must have unique non-empty reaction IDs.", call. = FALSE)
    }
    missing <- setdiff(reactions, penalty_names)
    unknown <- setdiff(penalty_names, reactions)
    if (length(missing)) {
      stop("Named `penalties` is missing reactions: ",
           paste(utils::head(missing, 10L), collapse = ", "), call. = FALSE)
    }
    if (length(unknown)) {
      stop("Named `penalties` contains unknown reactions: ",
           paste(utils::head(unknown, 10L), collapse = ", "), call. = FALSE)
    }
    penalty <- as.numeric(penalties[reactions])
  }
  if (any(!is.finite(penalty)) || any(penalty < 0)) {
    stop("`penalties` must contain finite non-negative values.", call. = FALSE)
  }
  penalty <- pmax(penalty, penalty_floor)

  identity <- Matrix::Diagonal(n)
  zero <- Matrix::Matrix(0, nrow = nrow(S), ncol = n, sparse = TRUE)
  mass_balance <- cbind(S, zero)
  abs_positive <- cbind(identity, -identity)
  abs_negative <- cbind(-identity, -identity)
  target <- Matrix::Matrix(0, nrow = 1L, ncol = 2L * n, sparse = TRUE)
  target[1L, target_index] <- if (identical(target_direction, "reverse")) -1 else 1

  absolute_upper <- pmax(abs(lb), abs(ub))
  if (any(!is.finite(absolute_upper))) {
    stop("Absolute flux bounds must be finite for the penalty LP.", call. = FALSE)
  }
  list(
    obj = c(rep(0, n), penalty),
    A = rbind(mass_balance, abs_positive, abs_negative, target),
    lhs = c(rep(0, nrow(S)), rep(-Inf, 2L * n), target_min),
    rhs = c(rep(0, nrow(S)), rep(0, 2L * n), Inf),
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

.rc_validate_lp_problem <- function(obj, A, lhs, rhs, lb, ub, time_limit) {
  if (is.null(dim(A)) || length(dim(A)) != 2L) {
    stop("`A` must be a two-dimensional constraint matrix.", call. = FALSE)
  }
  n_row <- nrow(A)
  n_col <- ncol(A)
  if (length(obj) != n_col || length(lb) != n_col || length(ub) != n_col) {
    stop("LP objective and variable bounds must align with ncol(A).", call. = FALSE)
  }
  if (length(lhs) != n_row || length(rhs) != n_row) {
    stop("LP row bounds must align with nrow(A).", call. = FALSE)
  }
  if (anyNA(obj) || any(!is.finite(obj))) {
    stop("LP objective coefficients must be finite.", call. = FALSE)
  }
  if (anyNA(lhs) || anyNA(rhs) || any(lhs > rhs)) {
    stop("LP row bounds must be non-missing and satisfy lhs <= rhs.", call. = FALSE)
  }
  if (anyNA(lb) || anyNA(ub) || any(!is.finite(lb)) || any(!is.finite(ub)) ||
      any(lb > ub)) {
    stop("LP variable bounds must be finite and satisfy lb <= ub.", call. = FALSE)
  }
  if (!is.numeric(time_limit) || length(time_limit) != 1L ||
      !is.finite(time_limit) || time_limit <= 0) {
    stop("`time_limit` must be one positive finite number.", call. = FALSE)
  }
  invisible(TRUE)
}

.rc_expand_ranged_constraints <- function(A, lhs, rhs) {
  A <- .rc_as_dgCMatrix(A)
  rows <- list()
  sense <- character()
  bound <- numeric()
  for (i in seq_len(nrow(A))) {
    lower_finite <- is.finite(lhs[[i]])
    upper_finite <- is.finite(rhs[[i]])
    if (!lower_finite && !upper_finite) next
    if (lower_finite && upper_finite && abs(lhs[[i]] - rhs[[i]]) <= 1e-12) {
      rows[[length(rows) + 1L]] <- A[i, , drop = FALSE]
      sense <- c(sense, "=")
      bound <- c(bound, lhs[[i]])
    } else {
      if (lower_finite) {
        rows[[length(rows) + 1L]] <- A[i, , drop = FALSE]
        sense <- c(sense, ">")
        bound <- c(bound, lhs[[i]])
      }
      if (upper_finite) {
        rows[[length(rows) + 1L]] <- A[i, , drop = FALSE]
        sense <- c(sense, "<")
        bound <- c(bound, rhs[[i]])
      }
    }
  }
  matrix_out <- if (length(rows)) do.call(rbind, rows) else
    Matrix::Matrix(0, nrow = 0L, ncol = ncol(A), sparse = TRUE)
  list(A = matrix_out, sense = sense, bound = bound)
}

rc_solve_lp <- function(obj, A, lhs, rhs, lb, ub,
                        solver = "highs", time_limit = 60) {
  solver <- match.arg(solver, c("highs", "gurobi", "glpk"))
  .rc_validate_lp_problem(obj, A, lhs, rhs, lb, ub, time_limit)
  if (solver == "gurobi") {
    return(rc_solve_lp_gurobi(obj, A, lhs, rhs, lb, ub, time_limit))
  }
  if (solver == "highs") {
    return(rc_solve_lp_highs(obj, A, lhs, rhs, lb, ub, time_limit))
  }
  rc_solve_lp_glpk(obj, A, lhs, rhs, lb, ub, time_limit)
}

rc_standardize_lp_result <- function(status, objective_value,
                                     primal_solution, runtime,
                                     message = NULL) {
  status_text <- tolower(trimws(as.character(status %||% "error")[[1L]]))
  standardized <- if (grepl("infeas", status_text)) {
    "infeasible"
  } else if (grepl("unbound", status_text)) {
    "unbounded"
  } else if (grepl("time|limit", status_text)) {
    "time_limit"
  } else if (!grepl("suboptimal|not[ _-]*optimal", status_text) &&
             grepl("(^|[^[:alpha:]])optimal([^[:alpha:]]|$)", status_text)) {
    "optimal"
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
  expanded <- .rc_expand_ranged_constraints(A, lhs, rhs)
  model <- list(
    modelsense = "min",
    obj = as.numeric(obj),
    A = expanded$A,
    sense = expanded$sense,
    rhs = expanded$bound,
    lb = as.numeric(lb),
    ub = as.numeric(ub)
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
  expanded <- .rc_expand_ranged_constraints(A, lhs, rhs)
  direction <- ifelse(expanded$sense == "=", "==",
                      ifelse(expanded$sense == ">", ">=", "<="))
  result <- tryCatch(
    Rglpk::Rglpk_solve_LP(
      obj, expanded$A, direction, expanded$bound,
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
  status <- if (result$status == 0) "optimal" else paste0("glpk_status_", result$status)
  rc_standardize_lp_result(
    status, result$optimum, result$solution,
    proc.time()[["elapsed"]] - started
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x
