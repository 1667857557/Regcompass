# LP representation and persistent solver runtime for exact Python CORDA2 semantics.

.rc_corda2_solver_feasibility_tolerance <- function(solver) {
  solver <- match.arg(as.character(solver), c("highs", "gurobi", "glpk"))
  switch(
    solver,
    highs = 1e-7,
    glpk = 1e-7,
    gurobi = 1e-6
  )
}

.rc_corda_split_model <- function(
    gem, tolerance = 1e-7, upper_bound = 1e6) {
  validated <- rc_validate_gem(gem)
  reactions <- validated$reactions
  reaction_lb <- as.numeric(validated$lb)
  reaction_ub <- as.numeric(validated$ub)
  names(reaction_lb) <- names(reaction_ub) <- reactions

  # Python CORDA.__init__ calls the COBRA bound setters in this exact order.
  # Each setter validates the transient lower/upper pair before the next setter.
  normalized_lb <- reaction_lb
  normalized_ub <- reaction_ub
  for (i in seq_along(reactions)) {
    reaction <- reactions[[i]]
    lower <- reaction_lb[[i]]
    upper <- reaction_ub[[i]]
    if (lower < -tolerance) {
      proposed_lower <- -upper_bound
      if (proposed_lower > upper) {
        stop(
          "Exact Python CORDA2 constructor would reject the transient lower ",
          "bound assignment for `", reaction, "`: ", proposed_lower,
          " > current upper bound ", upper, ".",
          call. = FALSE
        )
      }
      lower <- proposed_lower
    }
    if (upper > tolerance) {
      proposed_upper <- upper_bound
      if (lower > proposed_upper) {
        stop(
          "Exact Python CORDA2 constructor would reject the transient upper ",
          "bound assignment for `", reaction, "`: current lower bound ",
          lower, " > ", proposed_upper, ".",
          call. = FALSE
        )
      }
      upper <- proposed_upper
    }
    normalized_lb[[i]] <- lower
    normalized_ub[[i]] <- upper
  }

  n_reactions <- length(reactions)
  n_variables <- 2L * n_reactions
  columns <- vector("list", n_variables)
  direction_rows <- vector("list", n_variables)
  lower <- upper <- numeric(n_variables)
  variable_id <- character(n_variables)

  cursor <- 0L
  for (i in seq_along(reactions)) {
    reaction <- reactions[[i]]
    lb <- normalized_lb[[i]]
    ub <- normalized_ub[[i]]

    cursor <- cursor + 1L
    forward <- paste0(reaction, "::forward")
    variable_id[[cursor]] <- forward
    columns[[cursor]] <- validated$S[, i, drop = FALSE]
    if (lb > 0) {
      lower[[cursor]] <- lb
      upper[[cursor]] <- ub
    } else if (ub < 0) {
      lower[[cursor]] <- 0
      upper[[cursor]] <- 0
    } else {
      lower[[cursor]] <- 0
      upper[[cursor]] <- ub
    }
    direction_rows[[cursor]] <- data.frame(
      variable_id = forward,
      reaction_id = reaction,
      direction = "forward",
      original_index = i,
      stringsAsFactors = FALSE
    )

    cursor <- cursor + 1L
    reverse <- paste0(reaction, "::reverse")
    variable_id[[cursor]] <- reverse
    columns[[cursor]] <- -validated$S[, i, drop = FALSE]
    if (lb > 0) {
      lower[[cursor]] <- 0
      upper[[cursor]] <- 0
    } else if (ub < 0) {
      lower[[cursor]] <- -ub
      upper[[cursor]] <- -lb
    } else {
      lower[[cursor]] <- 0
      upper[[cursor]] <- -lb
    }
    direction_rows[[cursor]] <- data.frame(
      variable_id = reverse,
      reaction_id = reaction,
      direction = "reverse",
      original_index = i,
      stringsAsFactors = FALSE
    )
  }

  S <- do.call(cbind, columns)
  colnames(S) <- variable_id
  lower <- stats::setNames(as.numeric(lower), variable_id)
  upper <- stats::setNames(as.numeric(upper), variable_id)
  direction_table <- do.call(rbind, direction_rows)
  rownames(direction_table) <- NULL

  list(
    S = .rc_as_dgCMatrix(S),
    lb = lower,
    ub = upper,
    direction_table = direction_table,
    reaction_order = reactions,
    variable_to_reaction = stats::setNames(
      direction_table$reaction_id, direction_table$variable_id
    ),
    variable_to_direction = stats::setNames(
      direction_table$direction, direction_table$variable_id
    ),
    original_reaction_lb = reaction_lb,
    original_reaction_ub = reaction_ub,
    normalized_reaction_lb = normalized_lb,
    normalized_reaction_ub = normalized_ub,
    algorithm = "cobra_forward_reverse_variables_exact_CORDA2_bounds",
    tolerance = as.numeric(tolerance),
    upper_bound = as.numeric(upper_bound)
  )
}

.rc_corda_highs_api_available <- function() {
  if (!requireNamespace("highs", quietly = TRUE)) return(FALSE)
  required <- c(
    "highs_model", "hi_new_solver", "hi_solver_set_objective",
    "hi_solver_set_variable_bounds", "hi_solver_run",
    "hi_solver_status_message", "hi_solver_get_solution",
    "hi_solver_info", "hi_solver_set_option", "hi_solver_get_option"
  )
  all(required %in% getNamespaceExports("highs"))
}

.rc_corda_highs_call <- function(name, ...) {
  do.call(getExportedValue("highs", name), list(...))
}

.rc_corda_highs_exact_control <- function(tolerance, time_limit = Inf) {
  highs::highs_control(
    threads = 1L,
    time_limit = as.numeric(time_limit),
    log_to_console = FALSE,
    output_flag = FALSE,
    solver = "simplex",
    primal_feasibility_tolerance = as.numeric(tolerance)
  )
}

.rc_corda_highs_option_equal <- function(observed, expected) {
  if (is.numeric(expected)) {
    return(isTRUE(all.equal(
      as.numeric(observed), as.numeric(expected),
      tolerance = 0, check.attributes = FALSE
    )))
  }
  identical(as.character(observed), as.character(expected)) ||
    identical(as.logical(observed), as.logical(expected))
}

.rc_corda_highs_set_and_verify <- function(pointer, key, value) {
  .rc_corda_highs_call("hi_solver_set_option", pointer, key, value)
  observed <- .rc_corda_highs_call("hi_solver_get_option", pointer, key)
  if (!.rc_corda_highs_option_equal(observed, value)) {
    stop(
      "HiGHS did not retain the exact CORDA2 solver option `", key,
      "`: requested ", paste(value, collapse = ","),
      ", observed ", paste(observed, collapse = ","), ".",
      call. = FALSE
    )
  }
  observed
}

.rc_corda_changed_indices <- function(current, requested, label) {
  current <- as.numeric(current)
  requested <- as.numeric(requested)
  if (length(current) != length(requested)) {
    stop("CORDA2 ", label, " length changed during persistent solving.",
         call. = FALSE)
  }
  which(
    (current != requested) |
      (is.na(current) & !is.na(requested)) |
      (!is.na(current) & is.na(requested))
  )
}

.rc_corda_release_lp_engine <- function(engine) {
  if (!is.list(engine) || isTRUE(engine$released)) return(engine)
  if (!is.null(engine$pointer) && requireNamespace("highs", quietly = TRUE)) {
    exports <- getNamespaceExports("highs")
    if ("hi_solver_clear" %in% exports) {
      try(.rc_corda_highs_call("hi_solver_clear", engine$pointer), silent = TRUE)
    } else if ("hi_solver_clear_model" %in% exports) {
      try(.rc_corda_highs_call(
        "hi_solver_clear_model", engine$pointer
      ), silent = TRUE)
    }
  }
  engine$pointer <- NULL
  engine$released <- TRUE
  engine
}

.rc_corda_execution_metrics <- function(engine) {
  n_variables <- if (is.list(engine$split) && !is.null(engine$split$S)) {
    ncol(engine$split$S)
  } else {
    0L
  }
  full_values <- as.numeric(engine$n_full_vector_numeric_values %||% 0)
  transmitted <- as.numeric(engine$n_transmitted_numeric_values %||% 0)
  list(
    n_variables = as.integer(n_variables),
    n_solves = as.integer(engine$n_solves %||% 0L),
    n_fallback = as.integer(engine$n_fallback %||% 0L),
    n_objective_coeff_updates = as.integer(
      engine$n_objective_coeff_updates %||% 0L
    ),
    n_bound_index_updates = as.integer(
      engine$n_bound_index_updates %||% 0L
    ),
    n_sparse_update_calls = as.integer(
      engine$n_sparse_update_calls %||% 0L
    ),
    n_full_vector_numeric_values = full_values,
    n_transmitted_numeric_values = transmitted,
    n_full_vector_numeric_values_avoided = as.numeric(
      engine$n_full_vector_numeric_values_avoided %||% 0
    ),
    transmitted_fraction_of_full = if (full_values > 0) {
      transmitted / full_values
    } else {
      NA_real_
    },
    persistent_solver = identical(engine$type, "highs_persistent_cpp"),
    persistent_disabled = isTRUE(engine$persistent_disabled),
    solver_configuration_verified =
      isTRUE(engine$solver_configuration_verified),
    release_policy = "explicit_native_clear_on_reconstruction_exit"
  )
}

.rc_corda_new_lp_engine <- function(split, solver, time_limit) {
  n_variables <- ncol(split$S)
  engine <- list(
    type = "one_shot",
    pointer = NULL,
    split = split,
    solver = solver,
    time_limit = time_limit,
    n_solves = 0L,
    n_fallback = 0L,
    current_objective = rep(0, n_variables),
    current_lower = as.numeric(split$lb),
    current_upper = as.numeric(split$ub),
    n_objective_coeff_updates = 0L,
    n_bound_index_updates = 0L,
    n_sparse_update_calls = 0L,
    n_full_vector_numeric_values = 0,
    n_transmitted_numeric_values = 0,
    n_full_vector_numeric_values_avoided = 0,
    persistent_disabled = FALSE,
    solver_configuration_verified = FALSE,
    verified_solver_options = list(),
    released = FALSE
  )
  if (!identical(solver, "highs") || !.rc_corda_highs_api_available()) {
    return(engine)
  }

  persistent <- tryCatch({
    model <- .rc_corda_highs_call(
      "highs_model",
      L = rep(0, n_variables),
      lower = as.numeric(split$lb),
      upper = as.numeric(split$ub),
      A = split$S,
      lhs = rep(0, nrow(split$S)),
      rhs = rep(0, nrow(split$S)),
      maximum = FALSE
    )
    pointer <- .rc_corda_highs_call("hi_new_solver", model)
    settings <- list(
      output_flag = FALSE,
      threads = 1L,
      solver = "simplex",
      primal_feasibility_tolerance = as.numeric(split$tolerance)
    )
    if (is.finite(time_limit)) settings$time_limit <- as.numeric(time_limit)
    observed <- lapply(names(settings), function(key) {
      .rc_corda_highs_set_and_verify(pointer, key, settings[[key]])
    })
    names(observed) <- names(settings)
    list(pointer = pointer, observed = observed)
  }, error = function(e) e)

  if (inherits(persistent, "error")) {
    engine$persistent_disabled <- TRUE
    engine$solver_configuration_message <- conditionMessage(persistent)
    return(engine)
  }
  engine$type <- "highs_persistent_cpp"
  engine$pointer <- persistent$pointer
  engine$solver_configuration_verified <- TRUE
  engine$solver_configuration_message <-
    "verified by HiGHS option round trip"
  engine$verified_solver_options <- persistent$observed
  engine
}

.rc_corda_one_shot_solve <- function(
    engine, objective, lower, upper) {
  if (identical(engine$solver, "highs")) {
    answer <- tryCatch(
      highs::highs_solve(
        L = as.numeric(objective),
        lower = as.numeric(lower),
        upper = as.numeric(upper),
        A = engine$split$S,
        lhs = rep(0, nrow(engine$split$S)),
        rhs = rep(0, nrow(engine$split$S)),
        maximum = FALSE,
        control = .rc_corda_highs_exact_control(
          tolerance = engine$split$tolerance,
          time_limit = engine$time_limit
        )
      ),
      error = function(e) e
    )
    if (inherits(answer, "error")) {
      return(list(
        status = "error",
        solution = numeric(),
        objective = NA_real_,
        backend = "one_shot_highs_exact_configuration",
        solver_message = conditionMessage(answer)
      ))
    }
    return(list(
      status = .rc_lp_status(
        answer$status_message %||% answer$solver_msg %||% "",
        answer$status %||% NA_integer_
      ),
      solution = as.numeric(answer$primal_solution %||% numeric()),
      objective = as.numeric(answer$objective_value %||% NA_real_),
      backend = "one_shot_highs_exact_configuration",
      solver_message = answer$status_message %||% ""
    ))
  }

  answer <- rc_solve_lp(
    obj = objective,
    A = engine$split$S,
    lhs = rep(0, nrow(engine$split$S)),
    rhs = rep(0, nrow(engine$split$S)),
    lb = lower,
    ub = upper,
    solver = engine$solver,
    time_limit = engine$time_limit
  )
  list(
    status = answer$status,
    solution = answer$solution,
    objective = answer$objective,
    backend = "one_shot",
    solver_message = answer$solver_message %||% ""
  )
}

.rc_corda_engine_solve <- function(
    engine, objective, lower, upper) {
  objective <- as.numeric(objective)
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)
  n_variables <- length(objective)
  if (length(lower) != n_variables || length(upper) != n_variables) {
    stop("CORDA2 objective and bound vectors must have equal length.",
         call. = FALSE)
  }
  engine$n_solves <- engine$n_solves + 1L
  engine$n_full_vector_numeric_values <-
    engine$n_full_vector_numeric_values + 3 * n_variables

  if (!identical(engine$type, "highs_persistent_cpp") ||
      isTRUE(engine$persistent_disabled)) {
    engine$n_transmitted_numeric_values <-
      engine$n_transmitted_numeric_values + 3 * n_variables
    return(list(
      answer = .rc_corda_one_shot_solve(engine, objective, lower, upper),
      engine = engine
    ))
  }

  objective_index <- .rc_corda_changed_indices(
    engine$current_objective, objective, "objective"
  )
  bound_index <- union(
    .rc_corda_changed_indices(engine$current_lower, lower, "lower bound"),
    .rc_corda_changed_indices(engine$current_upper, upper, "upper bound")
  )
  transmitted <- length(objective_index) + 2L * length(bound_index)

  persistent <- tryCatch({
    if (length(objective_index)) {
      .rc_corda_highs_call(
        "hi_solver_set_objective", engine$pointer,
        index = as.integer(objective_index - 1L),
        coeff = objective[objective_index]
      )
    }
    if (length(bound_index)) {
      .rc_corda_highs_call(
        "hi_solver_set_variable_bounds", engine$pointer,
        index = as.integer(bound_index - 1L),
        lower = lower[bound_index],
        upper = upper[bound_index]
      )
    }
    .rc_corda_highs_call("hi_solver_run", engine$pointer)
    status_message <- .rc_corda_highs_call(
      "hi_solver_status_message", engine$pointer
    )
    status <- .rc_lp_status(status_message)
    solution <- if (identical(status, "optimal")) {
      as.numeric(.rc_corda_highs_call(
        "hi_solver_get_solution", engine$pointer
      )$col_value)
    } else {
      numeric()
    }
    info <- tryCatch(
      .rc_corda_highs_call("hi_solver_info", engine$pointer),
      error = function(e) list()
    )
    list(
      status = status,
      solution = solution,
      objective = as.numeric(
        info$objective_function_value %||% NA_real_
      ),
      backend = "highs_persistent_cpp_sparse_delta_basis_reuse",
      solver_message = as.character(status_message)
    )
  }, error = function(e) e)

  if (!inherits(persistent, "error")) {
    engine$current_objective <- objective
    engine$current_lower <- lower
    engine$current_upper <- upper
    engine$n_objective_coeff_updates <-
      engine$n_objective_coeff_updates + length(objective_index)
    engine$n_bound_index_updates <-
      engine$n_bound_index_updates + length(bound_index)
    engine$n_sparse_update_calls <- engine$n_sparse_update_calls + 1L
    engine$n_transmitted_numeric_values <-
      engine$n_transmitted_numeric_values + transmitted
    engine$n_full_vector_numeric_values_avoided <-
      engine$n_full_vector_numeric_values_avoided +
      (3 * n_variables - transmitted)
    return(list(answer = persistent, engine = engine))
  }

  failure_message <- conditionMessage(persistent)
  engine$n_fallback <- engine$n_fallback + 1L
  engine <- .rc_corda_release_lp_engine(engine)
  engine$type <- "one_shot"
  engine$persistent_disabled <- TRUE
  engine$current_objective <- objective
  engine$current_lower <- lower
  engine$current_upper <- upper
  engine$n_transmitted_numeric_values <-
    engine$n_transmitted_numeric_values + 3 * n_variables
  fallback <- .rc_corda_one_shot_solve(engine, objective, lower, upper)
  fallback$backend <- "highs_persistent_failed_one_shot_fallback"
  fallback$solver_message <- paste(
    failure_message, fallback$solver_message
  )
  list(answer = fallback, engine = engine)
}

