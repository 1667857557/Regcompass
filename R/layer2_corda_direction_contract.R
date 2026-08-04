# Match the original CORDA directional target semantics.

.rc_corda_target_bounds <- function(split, target, epsilon = NULL) {
  target <- as.character(target)
  if (length(target) != 1L || is.na(target) ||
      !target %in% split$direction_table$variable_id) {
    stop("CORDA target direction is not present in the split model.",
         call. = FALSE)
  }
  target_row <- split$direction_table[
    split$direction_table$variable_id == target, , drop = FALSE
  ]
  reaction <- as.character(target_row$reaction_id[[1L]])
  opposite <- as.character(split$direction_table$variable_id[
    split$direction_table$reaction_id == reaction &
      split$direction_table$variable_id != target
  ])
  lower <- split$lb
  upper <- split$ub
  if (length(opposite)) {
    lower[opposite] <- 0
    upper[opposite] <- 0
  }
  target_index <- match(target, names(upper))
  if (!is.null(epsilon)) {
    epsilon <- as.numeric(epsilon)
    if (length(epsilon) != 1L || !is.finite(epsilon) || epsilon <= 0) {
      stop("CORDA target epsilon must be one positive finite number.",
           call. = FALSE)
    }
    lower[[target_index]] <- max(lower[[target_index]], epsilon)
  }
  list(
    lower = lower,
    upper = upper,
    target = target,
    target_reaction = reaction,
    opposite_variables = opposite,
    target_index = target_index,
    original_code_semantics = paste(
      "test one signed direction of the original reaction while the",
      "opposite split direction is unavailable"
    )
  )
}

.rc_corda_dependency_task <- function(
    engine, task, confidence, options) {
  split <- engine$split
  target <- as.character(task$variable_id[[1L]])
  stage <- as.character(task$stage[[1L]])
  replicate <- as.integer(task$replicate[[1L]])
  noise_namespace <- as.character(options$noise_namespace %||% "")
  base_cost <- .rc_corda_base_cost(
    split, confidence, stage, options$gamma
  )
  objective <- as.numeric(base_cost) + .rc_corda_noise(
    length(base_cost), options$seed,
    c(noise_namespace, stage, target, replicate), options$kappa
  )
  if (!target %in% names(split$ub) || split$ub[[target]] < options$epsilon) {
    return(list(
      result = list(
        task = task, status = "target_blocked", associated = character(),
        target_flux = 0, objective = NA_real_, backend = engine$type,
        opposite_direction_blocked = character(),
        noise_namespace = noise_namespace
      ),
      engine = engine
    ))
  }
  bounds <- .rc_corda_target_bounds(
    split, target, epsilon = options$epsilon
  )
  solved <- .rc_corda_engine_solve(
    engine, objective, bounds$lower, bounds$upper
  )
  engine <- solved$engine
  answer <- solved$answer
  flux <- answer$solution
  associated <- character()
  target_flux <- NA_real_
  if (identical(answer$status, "optimal") &&
      length(flux) == ncol(split$S)) {
    names(flux) <- colnames(split$S)
    target_flux <- flux[[target]]
    penalized <- names(base_cost)[base_cost > 0]
    active <- penalized[flux[penalized] > options$flux_tolerance]
    associated <- unique(as.character(
      split$variable_to_reaction[active]
    ))
    associated <- setdiff(
      associated, as.character(task$reaction_id[[1L]])
    )
  }
  list(
    result = list(
      task = task,
      status = answer$status,
      associated = sort(associated),
      target_flux = target_flux,
      objective = answer$objective,
      backend = answer$backend,
      solver_message = answer$solver_message,
      opposite_direction_blocked = bounds$opposite_variables,
      noise_namespace = noise_namespace
    ),
    engine = engine
  )
}

.rc_corda_feasibility_task <- function(engine, task, options) {
  split <- engine$split
  target <- as.character(task$variable_id[[1L]])
  if (!target %in% names(split$ub) || split$ub[[target]] < options$epsilon) {
    return(list(
      result = list(
        task = task, status = "target_blocked", associated = character(),
        target_flux = 0, objective = NA_real_, backend = engine$type,
        opposite_direction_blocked = character()
      ),
      engine = engine
    ))
  }
  bounds <- .rc_corda_target_bounds(split, target, epsilon = NULL)
  objective <- rep(0, ncol(split$S))
  objective[[bounds$target_index]] <- -1
  solved <- .rc_corda_engine_solve(
    engine, objective, bounds$lower, bounds$upper
  )
  engine <- solved$engine
  answer <- solved$answer
  target_flux <- NA_real_
  status <- answer$status
  if (identical(status, "optimal") &&
      length(answer$solution) == ncol(split$S)) {
    target_flux <- as.numeric(answer$solution[[bounds$target_index]])
    if (!is.finite(target_flux) || target_flux < options$epsilon) {
      status <- "blocked"
    }
  }
  list(
    result = list(
      task = task,
      status = status,
      associated = character(),
      target_flux = target_flux,
      objective = answer$objective,
      backend = answer$backend,
      solver_message = answer$solver_message,
      opposite_direction_blocked = bounds$opposite_variables
    ),
    engine = engine
  )
}
