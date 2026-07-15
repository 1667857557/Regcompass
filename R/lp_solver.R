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
      rhs = numeric()
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
    rhs = bounds
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
