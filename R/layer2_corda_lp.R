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

  # Exact CORDA.__init__ bound normalization:
  # lower bounds below -tol become -UPPER and upper bounds above tol become UPPER.
  normalized_lb <- reaction_lb
  normalized_ub <- reaction_ub
  normalized_lb[normalized_lb < -tolerance] <- -upper_bound
  normalized_ub[normalized_ub > tolerance] <- upper_bound

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

.rc_corda_block_reactions <- function(split, reactions) {
  reactions <- unique(as.character(reactions))
  blocked <- split$direction_table$variable_id[
    split$direction_table$reaction_id %in% reactions
  ]
  split$ub[blocked] <- 0
  split$lb[blocked] <- pmin(split$lb[blocked], split$ub[blocked])
  split
}

.rc_corda_highs_api_available <- function() {
  if (!requireNamespace("highs", quietly = TRUE)) return(FALSE)
  required <- c(
    "highs_model", "hi_new_solver", "hi_solver_set_objective",
    "hi_solver_set_variable_bounds", "hi_solver_run",
    "hi_solver_status_message", "hi_solver_get_solution",
    "hi_solver_info", "hi_solver_set_option"
  )
  all(required %in% getNamespaceExports("highs"))
}

.rc_corda_highs_call <- function(name, ...) {
  do.call(getExportedValue("highs", name), list(...))
}

.rc_corda_new_lp_engine <- function(split, solver, time_limit) {
  if (identical(solver, "highs") && .rc_corda_highs_api_available()) {
    persistent <- tryCatch({
      model <- .rc_corda_highs_call(
        "highs_model",
        L = rep(0, ncol(split$S)),
        lower = as.numeric(split$lb),
        upper = as.numeric(split$ub),
        A = split$S,
        lhs = rep(0, nrow(split$S)),
        rhs = rep(0, nrow(split$S)),
        maximum = FALSE
      )
      pointer <- .rc_corda_highs_call("hi_new_solver", model)
      try(.rc_corda_highs_call(
        "hi_solver_set_option", pointer, "output_flag", "off"
      ), silent = TRUE)
      try(.rc_corda_highs_call(
        "hi_solver_set_option", pointer, "threads", "1"
      ), silent = TRUE)
      try(.rc_corda_highs_call(
        "hi_solver_set_option", pointer, "solver", "simplex"
      ), silent = TRUE)
      try(.rc_corda_highs_call(
        "hi_solver_set_option", pointer,
        "primal_feasibility_tolerance", as.character(split$tolerance)
      ), silent = TRUE)
      if (is.finite(time_limit)) {
        try(.rc_corda_highs_call(
          "hi_solver_set_option", pointer, "time_limit",
          as.character(time_limit)
        ), silent = TRUE)
      }
      list(
        type = "highs_persistent_cpp",
        pointer = pointer,
        split = split,
        solver = solver,
        time_limit = time_limit,
        n_solves = 0L,
        n_fallback = 0L
      )
    }, error = function(e) NULL)
    if (!is.null(persistent)) return(persistent)
  }
  list(
    type = "one_shot",
    pointer = NULL,
    split = split,
    solver = solver,
    time_limit = time_limit,
    n_solves = 0L,
    n_fallback = 0L
  )
}

.rc_corda_one_shot_solve <- function(
    engine, objective, lower, upper) {
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
  engine$n_solves <- engine$n_solves + 1L
  if (!identical(engine$type, "highs_persistent_cpp")) {
    return(list(
      answer = .rc_corda_one_shot_solve(engine, objective, lower, upper),
      engine = engine
    ))
  }
  persistent <- tryCatch({
    index <- seq_along(objective) - 1L
    .rc_corda_highs_call(
      "hi_solver_set_objective", engine$pointer,
      index = index, coeff = as.numeric(objective)
    )
    .rc_corda_highs_call(
      "hi_solver_set_variable_bounds", engine$pointer,
      index = index,
      lower = as.numeric(lower),
      upper = as.numeric(upper)
    )
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
      backend = "highs_persistent_cpp_basis_reuse",
      solver_message = as.character(status_message)
    )
  }, error = function(e) e)
  if (!inherits(persistent, "error")) {
    return(list(answer = persistent, engine = engine))
  }
  failure_message <- conditionMessage(persistent)
  engine$n_fallback <- engine$n_fallback + 1L
  fallback <- .rc_corda_one_shot_solve(engine, objective, lower, upper)
  fallback$backend <- "highs_persistent_failed_one_shot_fallback"
  fallback$solver_message <- paste(
    failure_message, fallback$solver_message
  )
  list(answer = fallback, engine = engine)
}

.rc_corda_task_bpparam <- function() {
  .rc_layer2_task_bpparam()
}

.rc_corda_worker_count <- function(BPPARAM, n_tasks) {
  # Exact CORDA2 reconstruction is serial within one model. Parallelism is
  # permitted only across independent cell-type x medium model instances.
  1L
}
