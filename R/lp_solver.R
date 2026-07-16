# Unified LP backends and directional COMPASS primitives.

.rc_lp_status <- function(message = "", code = NA_integer_) {
  text <- tolower(paste(message, collapse = " "))
  if (grepl("infeasible.*unbounded|unbounded.*infeasible", text)) {
    return("infeasible_or_unbounded")
  }
  if (grepl("infeasible", text)) return("infeasible")
  if (grepl("unbounded", text)) return("unbounded")
  if (grepl("time|limit", text)) return("time_limit")
  if (grepl("optimal", text)) return("optimal")
  if (is.finite(code) && as.integer(code) == 0L) return("optimal")
  "error"
}

.rc_expand_ranged_constraints <- function(A, lhs, rhs,
                                          tolerance = 1e-10) {
  if (!nrow(A)) {
    return(list(
      A = A,
      direction = character(),
      rhs = numeric(),
      sense = character(),
      bound = numeric()
    ))
  }

  source_rows <- integer()
  direction <- character()
  bounds <- numeric()

  for (i in seq_len(nrow(A))) {
    lower <- lhs[[i]]
    upper <- rhs[[i]]
    lower_finite <- is.finite(lower)
    upper_finite <- is.finite(upper)

    if (lower_finite && upper_finite &&
        abs(lower - upper) <= tolerance) {
      source_rows <- c(source_rows, i)
      direction <- c(direction, "==")
      bounds <- c(bounds, (lower + upper) / 2)
    } else {
      if (lower_finite) {
        source_rows <- c(source_rows, i)
        direction <- c(direction, ">=")
        bounds <- c(bounds, lower)
      }
      if (upper_finite) {
        source_rows <- c(source_rows, i)
        direction <- c(direction, "<=")
        bounds <- c(bounds, upper)
      }
    }
  }

  list(
    A = A[source_rows, , drop = FALSE],
    direction = direction,
    rhs = bounds,
    sense = ifelse(
      direction == "==", "=",
      ifelse(direction == ">=", ">", "<")
    ),
    bound = bounds
  )
}

rc_solve_lp <- function(obj, A, lhs, rhs, lb, ub,
                        solver = c("highs", "gurobi", "glpk"),
                        time_limit = 60) {
  solver <- match.arg(solver)
  obj <- as.numeric(obj)
  lb <- as.numeric(lb)
  ub <- as.numeric(ub)
  lhs <- as.numeric(lhs)
  rhs <- as.numeric(rhs)

  if (is.null(dim(A)) || length(dim(A)) != 2L) {
    stop("`A` must be a two-dimensional constraint matrix.", call. = FALSE)
  }
  if (ncol(A) != length(obj) || length(lb) != length(obj) ||
      length(ub) != length(obj)) {
    stop("LP objective, bounds and constraint columns are misaligned.",
         call. = FALSE)
  }
  if (nrow(A) != length(lhs) || nrow(A) != length(rhs)) {
    stop("LP constraint rows, lhs and rhs are misaligned.", call. = FALSE)
  }
  if (anyNA(obj) || any(!is.finite(obj)) || anyNA(lb) || anyNA(ub) ||
      anyNA(lhs) || anyNA(rhs) || any(lb > ub)) {
    stop("LP coefficients and bounds are invalid.", call. = FALSE)
  }
  if (!is.numeric(time_limit) || length(time_limit) != 1L ||
      is.na(time_limit) || time_limit <= 0) {
    stop("`time_limit` must be one positive number.", call. = FALSE)
  }

  A <- .rc_as_dgCMatrix(A)
  failure <- function(message) {
    list(
      status = "error",
      solution = numeric(),
      objective = NA_real_,
      solver = solver,
      solver_message = as.character(message)
    )
  }

  if (identical(solver, "highs")) {
    if (!requireNamespace("highs", quietly = TRUE)) {
      return(failure("Package 'highs' is not installed."))
    }
    answer <- tryCatch(
      highs::highs_solve(
        L = obj,
        lower = lb,
        upper = ub,
        A = A,
        lhs = lhs,
        rhs = rhs,
        maximum = FALSE,
        control = highs::highs_control(
          time_limit = time_limit,
          log_to_console = FALSE
        )
      ),
      error = function(e) e
    )
    if (inherits(answer, "error")) return(failure(conditionMessage(answer)))

    status <- .rc_lp_status(
      answer$status_message %||% answer$solver_msg %||% "",
      answer$status %||% NA_integer_
    )
    solution <- as.numeric(answer$primal_solution %||% numeric())
    objective <- as.numeric(answer$objective_value %||% NA_real_)
    return(list(
      status = status,
      solution = solution,
      objective = objective,
      solver = solver,
      solver_message = answer$status_message %||% ""
    ))
  }

  expanded <- .rc_expand_ranged_constraints(A, lhs, rhs)

  if (identical(solver, "glpk")) {
    if (!requireNamespace("Rglpk", quietly = TRUE)) {
      return(failure("Package 'Rglpk' is not installed."))
    }
    matrix_glpk <- expanded$A
    if (requireNamespace("slam", quietly = TRUE)) {
      matrix_glpk <- slam::as.simple_triplet_matrix(matrix_glpk)
    }
    answer <- tryCatch(
      Rglpk::Rglpk_solve_LP(
        obj = obj,
        mat = matrix_glpk,
        dir = expanded$direction,
        rhs = expanded$rhs,
        bounds = list(
          lower = list(ind = seq_along(lb), val = lb),
          upper = list(ind = seq_along(ub), val = ub)
        ),
        types = rep("C", length(obj)),
        max = FALSE,
        control = list(
          presolve = TRUE,
          tm_limit = if (is.finite(time_limit)) {
            as.integer(min(time_limit * 1000, .Machine$integer.max))
          } else {
            0L
          },
          canonicalize_status = TRUE
        )
      ),
      error = function(e) e
    )
    if (inherits(answer, "error")) return(failure(conditionMessage(answer)))

    status <- if (identical(as.integer(answer$status), 0L)) {
      "optimal"
    } else {
      "error"
    }
    return(list(
      status = status,
      solution = as.numeric(answer$solution %||% numeric()),
      objective = as.numeric(answer$objval %||% NA_real_),
      solver = solver,
      solver_message = paste0("Rglpk status ", answer$status)
    ))
  }

  if (!requireNamespace("gurobi", quietly = TRUE)) {
    return(failure("Package 'gurobi' is not installed."))
  }
  model <- list(
    A = expanded$A,
    obj = obj,
    lb = lb,
    ub = ub,
    sense = ifelse(expanded$direction == "==", "=",
                   ifelse(expanded$direction == ">=", ">", "<")),
    rhs = expanded$rhs,
    modelsense = "min"
  )
  parameters <- list(OutputFlag = 0)
  if (is.finite(time_limit)) parameters$TimeLimit <- time_limit
  answer <- tryCatch(
    gurobi::gurobi(model, params = parameters),
    error = function(e) e
  )
  if (inherits(answer, "error")) return(failure(conditionMessage(answer)))

  status <- .rc_lp_status(answer$status %||% "")
  list(
    status = status,
    solution = as.numeric(answer$x %||% numeric()),
    objective = as.numeric(answer$objval %||% NA_real_),
    solver = solver,
    solver_message = answer$status %||% ""
  )
}

rc_compass_vmax_directional <- function(S, lb, ub, target_reaction,
                                        direction = c("forward", "reverse"),
                                        solver = c("highs", "gurobi", "glpk"),
                                        time_limit = 60,
                                        flux_threshold = 1e-8) {
  direction <- match.arg(direction)
  solver <- match.arg(solver)
  reactions <- colnames(S)
  if (is.null(reactions) || !target_reaction %in% reactions) {
    stop("`target_reaction` is missing from the stoichiometric matrix.",
         call. = FALSE)
  }
  lb <- as.numeric(lb[reactions])
  ub <- as.numeric(ub[reactions])
  names(lb) <- names(ub) <- reactions
  target_index <- match(target_reaction, reactions)

  allowed <- if (identical(direction, "forward")) {
    ub[[target_index]] > flux_threshold
  } else {
    lb[[target_index]] < -flux_threshold
  }
  if (!isTRUE(allowed)) {
    return(list(
      feasible = FALSE,
      vmax = 0,
      status = "no_allowed_direction",
      flux = numeric()
    ))
  }

  objective <- rep(0, length(reactions))
  objective[[target_index]] <- if (identical(direction, "forward")) -1 else 1
  answer <- rc_solve_lp(
    obj = objective,
    A = S,
    lhs = rep(0, nrow(S)),
    rhs = rep(0, nrow(S)),
    lb = lb,
    ub = ub,
    solver = solver,
    time_limit = time_limit
  )
  if (!identical(answer$status, "optimal") ||
      length(answer$solution) != length(reactions)) {
    return(list(
      feasible = FALSE,
      vmax = NA_real_,
      status = answer$status,
      flux = numeric()
    ))
  }

  flux <- answer$solution
  names(flux) <- reactions
  vmax <- if (identical(direction, "forward")) {
    flux[[target_reaction]]
  } else {
    -flux[[target_reaction]]
  }
  feasible <- is.finite(vmax) && vmax >= flux_threshold
  list(
    feasible = feasible,
    vmax = if (is.finite(vmax)) max(0, vmax) else NA_real_,
    status = if (feasible) "optimal" else "blocked",
    flux = flux
  )
}

rc_compass_two_step_lp_directional <- function(
    S, lb, ub, target_reaction, penalties,
    target_direction = c("forward", "reverse"),
    omega = 0.95,
    solver = c("highs", "gurobi", "glpk"),
    time_limit = 60,
    flux_threshold = 1e-8) {
  target_direction <- match.arg(target_direction)
  solver <- match.arg(solver)
  if (!is.numeric(omega) || length(omega) != 1L ||
      !is.finite(omega) || omega <= 0 || omega > 1) {
    stop("`omega` must be in (0, 1].", call. = FALSE)
  }

  reactions <- colnames(S)
  if (is.null(reactions) || !target_reaction %in% reactions) {
    stop("`target_reaction` is missing from the stoichiometric matrix.",
         call. = FALSE)
  }
  lb <- as.numeric(lb[reactions])
  ub <- as.numeric(ub[reactions])
  names(lb) <- names(ub) <- reactions

  if (!is.null(names(penalties))) {
    missing_penalties <- setdiff(reactions, names(penalties))
    if (length(missing_penalties)) {
      stop("Reaction penalties are missing for: ",
           paste(utils::head(missing_penalties, 10L), collapse = ", "),
           call. = FALSE)
    }
    penalties <- as.numeric(penalties[reactions])
  } else {
    penalties <- as.numeric(penalties)
  }
  if (length(penalties) != length(reactions) ||
      any(!is.finite(penalties)) || any(penalties < 0)) {
    stop("`penalties` must provide one finite non-negative value per reaction.",
         call. = FALSE)
  }

  step1 <- rc_compass_vmax_directional(
    S = S,
    lb = lb,
    ub = ub,
    target_reaction = target_reaction,
    direction = target_direction,
    solver = solver,
    time_limit = time_limit,
    flux_threshold = flux_threshold
  )
  if (!isTRUE(step1$feasible)) {
    return(list(
      feasible = FALSE,
      penalty = NA_real_,
      vmax = step1$vmax,
      solver_status = step1$status,
      step1_status = step1$status,
      step2_status = "not_run",
      flux = numeric()
    ))
  }

  S <- .rc_as_dgCMatrix(S)
  n_reactions <- ncol(S)
  zero <- Matrix::Matrix(
    0,
    nrow = nrow(S),
    ncol = n_reactions,
    sparse = TRUE
  )
  mass_balance <- cbind(S, zero)

  positive <- Matrix::Matrix(
    0,
    nrow = n_reactions,
    ncol = 2L * n_reactions,
    sparse = TRUE
  )
  negative <- positive
  positive[cbind(seq_len(n_reactions), seq_len(n_reactions))] <- 1
  positive[cbind(seq_len(n_reactions), n_reactions + seq_len(n_reactions))] <- -1
  negative[cbind(seq_len(n_reactions), seq_len(n_reactions))] <- -1
  negative[cbind(seq_len(n_reactions), n_reactions + seq_len(n_reactions))] <- -1

  target <- Matrix::Matrix(
    0,
    nrow = 1,
    ncol = 2L * n_reactions,
    sparse = TRUE
  )
  target_index <- match(target_reaction, reactions)
  target[1, target_index] <- if (identical(target_direction, "forward")) 1 else -1

  A <- rbind(mass_balance, positive, negative, target)
  lhs <- c(
    rep(0, nrow(S)),
    rep(-Inf, 2L * n_reactions),
    omega * step1$vmax
  )
  rhs <- c(
    rep(0, nrow(S)),
    rep(0, 2L * n_reactions),
    Inf
  )
  auxiliary_upper <- pmax(abs(lb), abs(ub))

  step2 <- rc_solve_lp(
    obj = c(rep(0, n_reactions), penalties),
    A = A,
    lhs = lhs,
    rhs = rhs,
    lb = c(lb, rep(0, n_reactions)),
    ub = c(ub, auxiliary_upper),
    solver = solver,
    time_limit = time_limit
  )
  if (!identical(step2$status, "optimal") ||
      length(step2$solution) != 2L * n_reactions) {
    return(list(
      feasible = FALSE,
      penalty = NA_real_,
      vmax = step1$vmax,
      solver_status = step2$status,
      step1_status = step1$status,
      step2_status = step2$status,
      flux = numeric()
    ))
  }

  flux <- step2$solution[seq_len(n_reactions)]
  names(flux) <- reactions
  list(
    feasible = TRUE,
    penalty = max(0, as.numeric(step2$objective)),
    vmax = step1$vmax,
    solver_status = step2$status,
    step1_status = step1$status,
    step2_status = step2$status,
    flux = flux
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
