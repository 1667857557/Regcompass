# Exact paper-level CORDA costs and persistent-solver fallback semantics.

.rc_corda_base_cost <- function(split, confidence, stage, gamma) {
  reaction_confidence <- as.character(
    confidence[split$direction_table$reaction_id]
  )
  cost <- rep(0, nrow(split$direction_table))
  if (identical(stage, "stage1_hc_dependencies")) {
    cost[grepl("^MC", reaction_confidence)] <- sqrt(gamma)
    cost[reaction_confidence == "NC"] <- gamma
  } else if (identical(stage, "stage2_mc_nc_support")) {
    cost[reaction_confidence == "NC"] <- gamma
  } else if (identical(stage, "stage3_re_ot_dependencies")) {
    cost[reaction_confidence == "OT"] <- gamma
  } else {
    stop("Unknown CORDA dependency stage: ", stage, call. = FALSE)
  }
  stats::setNames(cost, split$direction_table$variable_id)
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
