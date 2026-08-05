# Verified solver configuration for exact Python CORDA2 execution semantics.

.rc_corda_highs_exact_control <- function(tolerance, time_limit = Inf) {
  args <- list(
    threads = 1L,
    time_limit = as.numeric(time_limit),
    log_to_console = FALSE,
    output_flag = FALSE,
    solver = "simplex",
    primal_feasibility_tolerance = as.numeric(tolerance)
  )
  do.call(highs::highs_control, args)
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
  .rc_corda_highs_call(
    "hi_solver_set_option", pointer, key, value
  )
  observed <- .rc_corda_highs_call(
    "hi_solver_get_option", pointer, key
  )
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

.rc_corda_new_lp_engine_verified_base <- .rc_corda_new_lp_engine

.rc_corda_new_lp_engine <- function(split, solver, time_limit) {
  engine <- .rc_corda_new_lp_engine_verified_base(
    split = split,
    solver = solver,
    time_limit = time_limit
  )
  if (!identical(engine$type, "highs_persistent_cpp")) return(engine)

  required <- c("hi_solver_set_option", "hi_solver_get_option")
  exports <- if (requireNamespace("highs", quietly = TRUE)) {
    getNamespaceExports("highs")
  } else {
    character()
  }
  if (!all(required %in% exports)) {
    engine <- .rc_corda_release_lp_engine(engine)
    engine$type <- "one_shot"
    engine$persistent_disabled <- TRUE
    engine$solver_configuration_verified <- FALSE
    engine$solver_configuration_message <-
      "HiGHS option round-trip API is unavailable; using exact one-shot path."
    return(engine)
  }

  settings <- list(
    output_flag = FALSE,
    threads = 1L,
    solver = "simplex",
    primal_feasibility_tolerance = as.numeric(split$tolerance)
  )
  if (is.finite(time_limit)) {
    settings$time_limit <- as.numeric(time_limit)
  }
  verified <- tryCatch({
    observed <- lapply(names(settings), function(key) {
      .rc_corda_highs_set_and_verify(
        engine$pointer, key, settings[[key]]
      )
    })
    names(observed) <- names(settings)
    observed
  }, error = function(e) e)

  if (inherits(verified, "error")) {
    message_text <- conditionMessage(verified)
    engine <- .rc_corda_release_lp_engine(engine)
    engine$type <- "one_shot"
    engine$persistent_disabled <- TRUE
    engine$solver_configuration_verified <- FALSE
    engine$solver_configuration_message <- message_text
    return(engine)
  }

  engine$solver_configuration_verified <- TRUE
  engine$solver_configuration_message <- "verified by HiGHS option round trip"
  engine$verified_solver_options <- verified
  engine
}

.rc_corda_one_shot_solve_verified_base <- .rc_corda_one_shot_solve

.rc_corda_one_shot_solve <- function(
    engine, objective, lower, upper) {
  if (!identical(engine$solver, "highs")) {
    return(.rc_corda_one_shot_solve_verified_base(
      engine = engine,
      objective = objective,
      lower = lower,
      upper = upper
    ))
  }
  if (!requireNamespace("highs", quietly = TRUE)) {
    return(list(
      status = "error",
      solution = numeric(),
      objective = NA_real_,
      backend = "one_shot_highs_exact_configuration",
      solver_message = "Package 'highs' is not installed."
    ))
  }

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
  list(
    status = .rc_lp_status(
      answer$status_message %||% answer$solver_msg %||% "",
      answer$status %||% NA_integer_
    ),
    solution = as.numeric(answer$primal_solution %||% numeric()),
    objective = as.numeric(answer$objective_value %||% NA_real_),
    backend = "one_shot_highs_exact_configuration",
    solver_message = answer$status_message %||% ""
  )
}
