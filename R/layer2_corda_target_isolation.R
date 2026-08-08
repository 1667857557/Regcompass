# Target-local solver lifecycle for deterministic CORDA2 stage parallelism.
#
# CORDA2 mathematical stages are independent across directional targets until
# their original stage barriers. HiGHS simplex basis state, however, is mutable.
# Reusing one solver instance across different targets can therefore make the
# selected primal optimum depend on target/chunk ordering when alternate optima
# exist. These helpers keep the LP formulation unchanged while limiting solver
# state reuse to the repeated solves belonging to one directional target.

.rc_corda_target_metric_engine <- function(split, solver, time_limit) {
  solver <- match.arg(as.character(solver), c("highs", "gurobi", "glpk"))
  n_variables <- ncol(split$S)
  list(
    type = if (identical(solver, "highs")) {
      "highs_persistent_cpp"
    } else {
      "one_shot"
    },
    pointer = NULL,
    released = TRUE,
    split = split,
    solver = solver,
    time_limit = time_limit,
    current_objective = rep(0, n_variables),
    current_lower = as.numeric(split$lb),
    current_upper = as.numeric(split$ub),
    n_solves = 0L,
    n_fallback = 0L,
    n_objective_coeff_updates = 0L,
    n_bound_index_updates = 0L,
    n_sparse_update_calls = 0L,
    n_full_vector_numeric_values = 0,
    n_transmitted_numeric_values = 0,
    n_full_vector_numeric_values_avoided = 0,
    persistent_disabled = FALSE,
    solver_configuration_verified = !identical(solver, "highs") || TRUE,
    target_engine_count = 0L,
    target_engine_scope = "fresh_solver_engine_per_directional_target",
    within_target_solver_state_reuse = TRUE,
    metric_accumulator_only = TRUE
  )
}

.rc_corda_new_target_engine <- function(engine, split = engine$split) {
  if (!is.list(engine) || is.null(engine$solver) || is.null(engine$time_limit)) {
    stop("CORDA2 target isolation requires a valid aggregate solver engine.",
         call. = FALSE)
  }
  target_engine <- .rc_corda_new_lp_engine(
    split = split,
    solver = engine$solver,
    time_limit = engine$time_limit
  )
  target_engine$regcompass_target_isolated <- TRUE
  target_engine
}

.rc_corda_absorb_target_engine <- function(aggregate_engine, target_engine) {
  if (!is.list(aggregate_engine) || !is.list(target_engine)) {
    stop("CORDA2 target engine accounting received an invalid engine.",
         call. = FALSE)
  }
  numeric_fields <- c(
    "n_solves",
    "n_fallback",
    "n_objective_coeff_updates",
    "n_bound_index_updates",
    "n_sparse_update_calls",
    "n_full_vector_numeric_values",
    "n_transmitted_numeric_values",
    "n_full_vector_numeric_values_avoided"
  )
  for (field in numeric_fields) {
    aggregate_engine[[field]] <-
      as.numeric(aggregate_engine[[field]] %||% 0) +
      as.numeric(target_engine[[field]] %||% 0)
  }
  aggregate_engine$persistent_disabled <-
    isTRUE(aggregate_engine$persistent_disabled) ||
    isTRUE(target_engine$persistent_disabled)
  if (identical(aggregate_engine$solver, "highs")) {
    aggregate_engine$solver_configuration_verified <-
      isTRUE(aggregate_engine$solver_configuration_verified) &&
      isTRUE(target_engine$solver_configuration_verified)
  }
  aggregate_engine$target_engine_count <-
    as.integer(aggregate_engine$target_engine_count %||% 0L) + 1L
  aggregate_engine$target_engine_scope <-
    "fresh_solver_engine_per_directional_target"
  aggregate_engine$within_target_solver_state_reuse <- TRUE
  target_engine <- .rc_corda_release_lp_engine(target_engine)
  aggregate_engine
}
