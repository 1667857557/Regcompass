# Weighted support LPs for CORDA-like Layer 2 completion.

.rc_corda_lp10 <- function(
    S, lb, ub, active_core, penalized_reactions, penalty_weights,
    epsilon, solver, time_limit, scaling_factor = 1e5) {
  active_core <- intersect(unique(as.character(active_core)), colnames(S))
  penalized_reactions <- intersect(
    unique(as.character(penalized_reactions)), colnames(S)
  )
  penalty_weights <- as.numeric(penalty_weights[penalized_reactions])
  if (anyNA(penalty_weights) || any(!is.finite(penalty_weights)) ||
      any(penalty_weights < 0)) {
    stop("CORDA-like support penalties must be finite and non-negative.",
         call. = FALSE)
  }
  n_reactions <- ncol(S)
  n_penalized <- length(penalized_reactions)
  if (!length(active_core)) {
    return(list(
      status = "empty_core", flux = numeric(), new_support = character(),
      objective = 0
    ))
  }
  S <- .rc_as_dgCMatrix(S)
  scaled_lb <- as.numeric(lb) * scaling_factor
  scaled_ub <- as.numeric(ub) * scaling_factor
  names(scaled_lb) <- names(lb)
  names(scaled_ub) <- names(ub)
  scaled_epsilon <- epsilon * scaling_factor
  zero <- Matrix::Matrix(
    0, nrow = nrow(S), ncol = n_penalized, sparse = TRUE
  )
  blocks <- list(cbind(S, zero))
  lhs <- rep(0, nrow(S))
  rhs <- rep(0, nrow(S))
  if (n_penalized) {
    positive <- Matrix::Matrix(
      0, nrow = n_penalized, ncol = n_reactions + n_penalized,
      sparse = TRUE
    )
    negative <- positive
    index <- match(penalized_reactions, colnames(S))
    positive[cbind(seq_len(n_penalized), index)] <- 1
    positive[cbind(
      seq_len(n_penalized), n_reactions + seq_len(n_penalized)
    )] <- -1
    negative[cbind(seq_len(n_penalized), index)] <- -1
    negative[cbind(
      seq_len(n_penalized), n_reactions + seq_len(n_penalized)
    )] <- -1
    blocks <- c(blocks, list(positive, negative))
    lhs <- c(lhs, rep(-Inf, 2L * n_penalized))
    rhs <- c(rhs, rep(0, 2L * n_penalized))
  }
  core_constraint <- Matrix::Matrix(
    0, nrow = length(active_core),
    ncol = n_reactions + n_penalized, sparse = TRUE
  )
  core_constraint[cbind(
    seq_along(active_core), match(active_core, colnames(S))
  )] <- 1
  blocks <- c(blocks, list(core_constraint))
  lhs <- c(lhs, rep(scaled_epsilon, length(active_core)))
  rhs <- c(rhs, rep(Inf, length(active_core)))
  auxiliary_upper <- if (n_penalized) {
    pmax(
      abs(scaled_lb[penalized_reactions]),
      abs(scaled_ub[penalized_reactions])
    )
  } else {
    numeric()
  }
  answer <- rc_solve_lp(
    obj = c(rep(0, n_reactions), penalty_weights),
    A = do.call(rbind, blocks),
    lhs = lhs,
    rhs = rhs,
    lb = c(scaled_lb, rep(0, n_penalized)),
    ub = c(scaled_ub, auxiliary_upper),
    solver = solver,
    time_limit = time_limit
  )
  if (!identical(answer$status, "optimal")) {
    return(list(
      status = answer$status, flux = numeric(), new_support = character(),
      objective = NA_real_
    ))
  }
  scaled_flux <- answer$solution[seq_len(n_reactions)]
  flux <- scaled_flux / scaling_factor
  names(flux) <- colnames(S)
  scaled_auxiliary <- if (n_penalized) {
    answer$solution[n_reactions + seq_len(n_penalized)]
  } else {
    numeric()
  }
  support_tolerance_scaled <- max(1e-7, scaled_epsilon * 1e-9)
  new_support <- penalized_reactions[
    scaled_auxiliary > support_tolerance_scaled
  ]
  list(
    status = answer$status,
    flux = flux,
    new_support = new_support,
    objective = answer$objective / scaling_factor,
    scaling_factor = scaling_factor,
    support_tolerance = support_tolerance_scaled / scaling_factor
  )
}

.rc_corda_complete_direction <- function(
    parent, biological_reactions, selected_support, targets, direction,
    epsilon, solver, time_limit, max_support_reactions, support_costs,
    scaling_factor = 1e5) {
  if (!nrow(targets)) {
    return(list(
      support = selected_support,
      unresolved = targets,
      iterations = data.frame()
    ))
  }
  validated <- rc_validate_gem(parent)
  reverse_targets <- if (identical(direction, "reverse")) {
    as.character(targets$reaction_id)
  } else {
    character()
  }
  oriented <- .rc_orient_reactions(
    validated$S, validated$lb, validated$ub, reverse_targets
  )
  remaining <- unique(as.character(targets$reaction_id))
  iteration_rows <- list()
  iteration <- 0L
  local_feasible <- function(reactions, core_ids) {
    local <- .rc_subset_gem(
      list(S = oriented$S, lb = oriented$lb, ub = oriented$ub), reactions
    )
    do.call(rbind, lapply(core_ids, function(reaction) {
      answer <- rc_compass_vmax_directional(
        local$S, local$lb, local$ub,
        reaction,
        direction = "forward",
        solver = solver,
        time_limit = time_limit,
        flux_threshold = epsilon
      )
      data.frame(
        reaction_id = reaction,
        feasible = isTRUE(answer$feasible),
        stringsAsFactors = FALSE
      )
    }))
  }
  while (length(remaining)) {
    iteration <- iteration + 1L
    before <- remaining
    unpenalized <- union(biological_reactions, selected_support)
    penalized <- setdiff(colnames(oriented$S), unpenalized)
    lp7 <- .rc_fastcore_lp7(
      oriented$S, oriented$lb, oriented$ub,
      remaining, epsilon, solver, time_limit
    )
    active <- lp7$active
    singleton_mode <- FALSE
    lp10_status <- "not_run"
    objective <- NA_real_
    added <- character()
    if (length(active)) {
      lp10 <- .rc_corda_lp10(
        oriented$S, oriented$lb, oriented$ub,
        active, penalized, support_costs,
        epsilon, solver, time_limit,
        scaling_factor = scaling_factor
      )
      lp10_status <- lp10$status
      objective <- lp10$objective
      added <- lp10$new_support
    }
    if (!length(active) || !identical(lp10_status, "optimal")) {
      singleton_mode <- TRUE
      for (reaction in remaining) {
        penalized <- setdiff(
          colnames(oriented$S), union(biological_reactions, selected_support)
        )
        one <- .rc_corda_lp10(
          oriented$S, oriented$lb, oriented$ub,
          reaction, penalized, support_costs,
          epsilon, solver, time_limit,
          scaling_factor = scaling_factor
        )
        if (identical(one$status, "optimal")) {
          selected_support <- union(selected_support, one$new_support)
          added <- union(added, one$new_support)
        }
      }
    } else {
      selected_support <- union(selected_support, added)
    }
    if (length(selected_support) > max_support_reactions) {
      stop(
        "CORDA-like support completion exceeded `max_support_reactions`.",
        call. = FALSE
      )
    }
    current_set <- union(biological_reactions, selected_support)
    check <- local_feasible(current_set, remaining)
    remaining <- as.character(check$reaction_id[!check$feasible])
    repaired <- setdiff(before, remaining)
    iteration_rows[[iteration]] <- data.frame(
      iteration = iteration,
      target_direction = direction,
      n_targets_before = length(before),
      n_targets_lp7_active = length(active),
      n_targets_repaired = length(repaired),
      n_new_support = length(added),
      lp7_status = lp7$status,
      weighted_lp10_status = lp10_status,
      weighted_support_objective = objective,
      epsilon = epsilon,
      lp10_scaling_factor = scaling_factor,
      singleton_mode = singleton_mode,
      stringsAsFactors = FALSE
    )
    if (!length(repaired)) break
  }
  unresolved <- targets[
    as.character(targets$reaction_id) %in% remaining,
    , drop = FALSE
  ]
  list(
    support = selected_support,
    unresolved = unresolved,
    iterations = if (length(iteration_rows)) {
      do.call(rbind, iteration_rows)
    } else {
      data.frame()
    }
  )
}
