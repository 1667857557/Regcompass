# Exact-source audit contract.
# No mathematical correction is applied to the pinned Python CORDA2 behavior.

# Python CORDA.__init__ updates reaction bounds through two ordered COBRA
# property setters. Each setter validates the transient reaction bounds before
# the second assignment occurs. Validate that same transient program before the
# split-variable representation is constructed.
.rc_corda_split_model_exact_base <- .rc_corda_split_model

.rc_corda2_validate_constructor_bound_order <- function(
    gem, tolerance = 1e-7, upper_bound = 1e6) {
  validated <- rc_validate_gem(gem)
  for (i in seq_along(validated$reactions)) {
    reaction <- validated$reactions[[i]]
    lower <- as.numeric(validated$lb[[i]])
    upper <- as.numeric(validated$ub[[i]])

    # Source order:
    # if r.lower_bound < -tol: r.lower_bound = -UPPER
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

    # Source order:
    # if r.upper_bound > tol: r.upper_bound = UPPER
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
  }
  invisible(TRUE)
}

.rc_corda_split_model <- function(
    gem, tolerance = 1e-7, upper_bound = 1e6) {
  .rc_corda2_validate_constructor_bound_order(
    gem = gem,
    tolerance = tolerance,
    upper_bound = upper_bound
  )
  .rc_corda_split_model_exact_base(
    gem = gem,
    tolerance = tolerance,
    upper_bound = upper_bound
  )
}

.rc_corda_build_three_stage_exact_base <- .rc_corda_build_three_stage

.rc_corda_build_three_stage <- function(
    split, classes, options, solver, time_limit = Inf) {
  requested_time_limit <- time_limit
  # Python CORDA exposes no time-limit constructor control and calls the active
  # solver without adding one. A generic RegCompass stage timeout must not turn
  # a slow target into an algorithmic `impossible` target.
  answer <- .rc_corda_build_three_stage_exact_base(
    split = split,
    classes = classes,
    options = options,
    solver = solver,
    time_limit = Inf
  )
  answer$source_fidelity <- "exact_for_met_prod_NULL"
  answer$intentional_corrections <- character()
  answer$solver_time_limit <- Inf
  answer$requested_regcompass_time_limit_ignored <- requested_time_limit
  answer
}
