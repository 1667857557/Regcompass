# Exact execution semantics of resendislab/corda Python CORDA2.
# Reference: c02e06d50606bf93f23d8f2e6d6ade0e996ca70e.

.rc_corda2_directional_confidence <- function(split, classes) {
  reaction_confidence <- stats::setNames(
    rep(0L, length(classes$confidence)), names(classes$confidence)
  )
  reaction_confidence[classes$hc] <- 3L
  reaction_confidence[classes$mc_module] <- 2L
  reaction_confidence[classes$mc_evidence] <- 1L
  reaction_confidence[classes$nc] <- -1L
  reaction_confidence[classes$ot] <- 0L
  value <- as.integer(
    reaction_confidence[split$direction_table$reaction_id]
  )
  if (anyNA(value) || any(!value %in% c(-1L, 0L, 1L, 2L, 3L))) {
    stop("Every CORDA2 reaction requires confidence in {-1,0,1,2,3}.",
         call. = FALSE)
  }
  stats::setNames(value, split$direction_table$variable_id)
}

.rc_corda2_forward_variable <- function(split, reaction) {
  variable <- split$direction_table$variable_id[
    split$direction_table$reaction_id == reaction &
      split$direction_table$direction == "forward"
  ]
  if (length(variable) != 1L) {
    stop("CORDA2 split model has no unique forward variable for `",
         reaction, "`.", call. = FALSE)
  }
  variable
}

.rc_corda2_penalties <- function(
    split, directional_confidence, penalize_medium, penalty_factor) {
  penalty <- stats::setNames(
    rep(0, ncol(split$S)), colnames(split$S)
  )
  penalized <- character()
  for (reaction in split$reaction_order) {
    forward <- .rc_corda2_forward_variable(split, reaction)
    confidence <- directional_confidence[[forward]]
    matched <- FALSE
    value <- if (isTRUE(penalize_medium) && confidence %in% c(1L, 2L)) {
      matched <- TRUE
      1
    } else if (identical(confidence, -1L)) {
      matched <- TRUE
      penalty_factor
    } else {
      0
    }
    if (isTRUE(matched)) {
      variables <- split$direction_table$variable_id[
        split$direction_table$reaction_id == reaction
      ]
      penalty[variables] <- value
      penalized <- c(penalized, variables)
    }
  }
  attr(penalty, "penalized_variables") <- unique(penalized)
  penalty
}

.rc_corda2_target_result <- function(
    split, target, stage, kind, status, associated = character(),
    target_flux = NA_real_, objective = NA_real_, backend = "",
    solver_message = "", redundancies = 0L, n_solves = 0L) {
  row <- split$direction_table[
    split$direction_table$variable_id == target, , drop = FALSE
  ]
  list(
    target = target,
    reaction_id = as.character(row$reaction_id[[1L]]),
    direction = as.character(row$direction[[1L]]),
    stage = stage,
    kind = kind,
    status = status,
    associated = as.character(associated),
    target_flux = as.numeric(target_flux),
    objective = as.numeric(objective),
    backend = as.character(backend),
    solver_message = as.character(solver_message),
    opposite_direction_blocked = character(),
    redundancies = as.integer(redundancies),
    n_solves = as.integer(n_solves)
  )
}

.rc_corda2_associated <- function(
    engine, split, targets, directional_confidence, options,
    penalize_medium = TRUE, redundancies = TRUE, stage) {
  targets <- as.character(targets)
  if (anyNA(targets) || any(!targets %in% colnames(split$S))) {
    stop("CORDA2 target variable is missing from the solver model.",
         call. = FALSE)
  }
  penalties <- .rc_corda2_penalties(
    split = split,
    directional_confidence = directional_confidence,
    penalize_medium = penalize_medium,
    penalty_factor = options$penalty_factor
  )
  penalized_variables <- attr(penalties, "penalized_variables") %||%
    character()
  max_iter <- if (isTRUE(redundancies)) options$redundancies else 1L
  needed_all <- character()
  impossible <- character()
  results <- vector("list", length(targets))
  redundancy_map <- stats::setNames(integer(length(targets)), targets)

  for (i in seq_along(targets)) {
    target <- targets[[i]]
    if (split$ub[[target]] < split$tolerance) {
      directional_confidence[[target]] <- -1L
      impossible <- c(impossible, target)
      results[[i]] <- .rc_corda2_target_result(
        split, target, stage, "dependency", "target_blocked",
        target_flux = 0,
        backend = engine$type,
        solver_message = "target upper bound is below CORDA2 tolerance"
      )
      next
    }

    bounds <- .rc_corda_target_bounds(
      split, target, epsilon = options$target_flux
    )
    penalty <- penalties
    needed_for_target <- character()
    has_new <- TRUE
    iteration <- 0L
    redundancy_count <- 0L
    status <- "not_run"
    target_flux <- NA_real_
    objective_value <- NA_real_
    backend <- engine$type
    solver_message <- ""

    while (isTRUE(has_new) && iteration < max_iter) {
      solved <- .rc_corda_engine_solve(
        engine,
        objective = as.numeric(penalty),
        lower = bounds$lower,
        upper = bounds$upper
      )
      engine <- solved$engine
      answer <- solved$answer
      iteration <- iteration + 1L
      status <- answer$status
      objective_value <- answer$objective
      backend <- answer$backend
      solver_message <- answer$solver_message %||% ""
      if (!identical(status, "optimal") ||
          length(answer$solution) != ncol(split$S)) {
        directional_confidence[[target]] <- -1L
        impossible <- c(impossible, target)
        break
      }

      flux <- as.numeric(answer$solution)
      names(flux) <- colnames(split$S)
      target_flux <- flux[[target]]
      need <- names(flux)[
        flux > split$tolerance &
          directional_confidence[names(flux)] %in% c(-1L, 1L, 2L) &
          names(flux) != target
      ]
      new <- need[!need %in% needed_for_target]
      has_new <- length(new) > 0L
      if (isTRUE(redundancies)) {
        redundancy_count <- redundancy_count + as.integer(has_new)
      }
      weighted_new <- new[new %in% penalized_variables]
      if (length(weighted_new)) {
        penalty[weighted_new] <- penalty[weighted_new] * options$cost_increase
      }
      needed_for_target <- sort(unique(c(needed_for_target, need)), method = "radix")
    }

    redundancy_map[[target]] <- redundancy_count
    needed_all <- c(needed_all, needed_for_target)
    results[[i]] <- .rc_corda2_target_result(
      split, target, stage, "dependency", status,
      associated = needed_for_target,
      target_flux = target_flux,
      objective = objective_value,
      backend = backend,
      solver_message = solver_message,
      redundancies = redundancy_count,
      n_solves = iteration
    )
  }

  list(
    engine = engine,
    confidence = directional_confidence,
    needed = needed_all,
    results = results,
    impossible = impossible,
    redundancies = redundancy_map,
    execution = list(
      n_targets = length(targets),
      n_chunks = if (length(targets)) 1L else 0L,
      workers = 1L,
      task_granularity = "python_serial_target_order",
      stage_barrier = TRUE,
      target_parallelism = FALSE,
      persistent_solver = identical(engine$type, "highs_persistent_cpp"),
      solver_runtime = engine$type,
      n_solves = engine$n_solves,
      n_fallback = engine$n_fallback
    )
  )
}
