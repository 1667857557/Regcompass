# FASTCC and add-only FASTCORE reconstruction for GRN-defined meta-modules.

.rc_subset_gem <- function(gem, reactions) {
  validated <- rc_validate_gem(gem)
  keep <- intersect(unique(as.character(reactions)), validated$reactions)
  if (!length(keep)) {
    stop("No reactions remain in the requested GEM subset.", call. = FALSE)
  }

  output <- gem
  output$S <- validated$S[, keep, drop = FALSE]
  output$lb <- validated$lb[keep]
  output$ub <- validated$ub[keep]
  used_metabolites <- rownames(output$S)[
    Matrix::rowSums(abs(output$S) > 0) > 0
  ]
  output$S <- output$S[used_metabolites, , drop = FALSE]

  if (!is.null(gem$reaction_meta) && is.data.frame(gem$reaction_meta)) {
    output$reaction_meta <- gem$reaction_meta[
      match(keep, as.character(gem$reaction_meta$reaction_id)),
      , drop = FALSE
    ]
  }
  if (!is.null(gem$metabolite_meta) && is.data.frame(gem$metabolite_meta)) {
    output$metabolite_meta <- gem$metabolite_meta[
      match(
        used_metabolites,
        as.character(gem$metabolite_meta$metabolite_id)
      ),
      , drop = FALSE
    ]
  }
  if (!is.null(gem$gpr_table) && is.data.frame(gem$gpr_table)) {
    output$gpr_table <- gem$gpr_table[
      as.character(gem$gpr_table$reaction_id) %in% keep,
      , drop = FALSE
    ]
  }
  output
}

.rc_orient_reactions <- function(S, lb, ub, reverse_reactions = character()) {
  reverse_reactions <- intersect(
    unique(as.character(reverse_reactions)),
    colnames(S)
  )
  if (!length(reverse_reactions)) {
    return(list(S = S, lb = lb, ub = ub))
  }
  index <- match(reverse_reactions, colnames(S))
  S[, index] <- -S[, index, drop = FALSE]
  old_lb <- lb[index]
  old_ub <- ub[index]
  lb[index] <- -old_ub
  ub[index] <- -old_lb
  list(S = S, lb = lb, ub = ub)
}

.rc_directional_feasibility <- function(gem, targets, solver = "highs",
                                        time_limit = 60,
                                        flux_threshold = 1e-8) {
  required <- c("reaction_id", "target_direction")
  if (!is.data.frame(targets) || !all(required %in% colnames(targets))) {
    stop("`targets` must contain reaction_id and target_direction.", call. = FALSE)
  }
  if (!nrow(targets)) {
    return(data.frame(
      reaction_id = character(), target_direction = character(),
      feasible = logical(), vmax = numeric(), solver_status = character(),
      stringsAsFactors = FALSE
    ))
  }
  validated <- rc_validate_gem(gem)
  rows <- lapply(seq_len(nrow(targets)), function(i) {
    reaction <- as.character(targets$reaction_id[[i]])
    direction <- as.character(targets$target_direction[[i]])
    if (!reaction %in% validated$reactions) {
      return(data.frame(
        reaction_id = reaction, target_direction = direction,
        feasible = FALSE, vmax = NA_real_, solver_status = "reaction_missing",
        stringsAsFactors = FALSE
      ))
    }
    if (!direction %in% c("forward", "reverse")) {
      return(data.frame(
        reaction_id = reaction, target_direction = direction,
        feasible = FALSE, vmax = 0, solver_status = "no_allowed_direction",
        stringsAsFactors = FALSE
      ))
    }
    answer <- rc_compass_vmax_directional(
      S = validated$S, lb = validated$lb, ub = validated$ub,
      target_reaction = reaction, direction = direction,
      solver = solver, time_limit = time_limit,
      flux_threshold = flux_threshold
    )
    data.frame(
      reaction_id = reaction, target_direction = direction,
      feasible = isTRUE(answer$feasible), vmax = answer$vmax,
      solver_status = answer$status, stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# FASTCORE LP-7: maximize the number of core reactions carrying epsilon flux.
.rc_fastcore_lp7 <- function(S, lb, ub, core_reactions, epsilon,
                             solver, time_limit) {
  core_reactions <- intersect(
    unique(as.character(core_reactions)),
    colnames(S)
  )
  n_reactions <- ncol(S)
  n_core <- length(core_reactions)
  if (!n_core) {
    return(list(
      status = "empty_core", flux = numeric(),
      active = character(), objective = 0
    ))
  }

  S <- .rc_as_dgCMatrix(S)
  zero <- Matrix::Matrix(
    0, nrow = nrow(S), ncol = n_core, sparse = TRUE
  )
  mass_balance <- cbind(S, zero)
  activation <- Matrix::Matrix(
    0, nrow = n_core, ncol = n_reactions + n_core,
    sparse = TRUE
  )
  activation[cbind(
    seq_len(n_core),
    match(core_reactions, colnames(S))
  )] <- 1
  activation[cbind(
    seq_len(n_core),
    n_reactions + seq_len(n_core)
  )] <- -1

  answer <- rc_solve_lp(
    obj = c(rep(0, n_reactions), rep(-1, n_core)),
    A = rbind(mass_balance, activation),
    lhs = c(rep(0, nrow(S)), rep(0, n_core)),
    rhs = c(rep(0, nrow(S)), rep(Inf, n_core)),
    lb = c(lb, rep(0, n_core)),
    ub = c(ub, rep(epsilon, n_core)),
    solver = solver,
    time_limit = time_limit
  )
  if (!identical(answer$status, "optimal")) {
    return(list(
      status = answer$status, flux = numeric(),
      active = character(), objective = NA_real_
    ))
  }
  flux <- answer$solution[seq_len(n_reactions)]
  names(flux) <- colnames(S)
  z <- answer$solution[n_reactions + seq_len(n_core)]
  active <- core_reactions[z >= epsilon * (1 - 1e-7)]
  list(
    status = answer$status,
    flux = flux,
    active = active,
    objective = -answer$objective
  )
}

.rc_fastcc_consistent_reactions <- function(gem, solver = "highs",
                                            time_limit = 300,
                                            epsilon = 1e-4) {
  validated <- rc_validate_gem(gem)
  if (!is.finite(epsilon) || epsilon <= 0) {
    stop("`epsilon` must be positive.", call. = FALSE)
  }

  collect_active <- function(S, lb, ub, candidates) {
    candidates <- unique(as.character(candidates))
    active_all <- character()
    remaining <- candidates
    while (length(remaining)) {
      batch <- .rc_fastcore_lp7(
        S, lb, ub, remaining, epsilon,
        solver = solver, time_limit = time_limit
      )
      active <- batch$active
      if (!length(active)) break
      active_all <- union(active_all, active)
      remaining <- setdiff(remaining, active)
    }
    if (length(remaining)) {
      for (reaction in remaining) {
        one <- .rc_fastcore_lp7(
          S, lb, ub, reaction, epsilon,
          solver = solver, time_limit = time_limit
        )
        if (length(one$active)) {
          active_all <- union(active_all, reaction)
        }
      }
    }
    active_all
  }

  forward_candidates <- validated$reactions[
    validated$ub >= epsilon
  ]
  consistent <- collect_active(
    validated$S, validated$lb, validated$ub,
    forward_candidates
  )

  reverse_candidates <- setdiff(
    validated$reactions[validated$lb <= -epsilon],
    consistent
  )
  if (length(reverse_candidates)) {
    oriented <- .rc_orient_reactions(
      validated$S, validated$lb, validated$ub,
      reverse_candidates
    )
    consistent <- union(
      consistent,
      collect_active(
        oriented$S, oriented$lb, oriented$ub,
        reverse_candidates
      )
    )
  }
  consistent
}

# FASTCORE LP-10. Following Vlassis et al., the core constraint and all flux
# bounds are scaled by 1e5; support is extracted using the original epsilon.
.rc_fastcore_lp10 <- function(S, lb, ub, active_core,
                              penalized_reactions,
                              epsilon, solver, time_limit,
                              scaling_factor = 1e5) {
  active_core <- intersect(
    unique(as.character(active_core)),
    colnames(S)
  )
  penalized_reactions <- intersect(
    unique(as.character(penalized_reactions)),
    colnames(S)
  )
  n_reactions <- ncol(S)
  n_penalized <- length(penalized_reactions)
  if (!length(active_core)) {
    return(list(
      status = "empty_core", flux = numeric(),
      new_support = character(), objective = 0
    ))
  }
  if (!is.finite(epsilon) || epsilon <= 0) {
    stop("`epsilon` must be positive.", call. = FALSE)
  }
  if (!is.finite(scaling_factor) || scaling_factor < 1) {
    stop("`scaling_factor` must be >= 1.", call. = FALSE)
  }

  S <- .rc_as_dgCMatrix(S)
  scaled_lb <- as.numeric(lb) * scaling_factor
  scaled_ub <- as.numeric(ub) * scaling_factor
  names(scaled_lb) <- names(lb)
  names(scaled_ub) <- names(ub)
  scaled_epsilon <- epsilon * scaling_factor

  zero <- Matrix::Matrix(
    0, nrow = nrow(S), ncol = n_penalized,
    sparse = TRUE
  )
  blocks <- list(cbind(S, zero))
  lhs <- rep(0, nrow(S))
  rhs <- rep(0, nrow(S))

  if (n_penalized) {
    positive <- Matrix::Matrix(
      0,
      nrow = n_penalized,
      ncol = n_reactions + n_penalized,
      sparse = TRUE
    )
    negative <- positive
    penalized_index <- match(penalized_reactions, colnames(S))
    positive[cbind(seq_len(n_penalized), penalized_index)] <- 1
    positive[cbind(
      seq_len(n_penalized),
      n_reactions + seq_len(n_penalized)
    )] <- -1
    negative[cbind(seq_len(n_penalized), penalized_index)] <- -1
    negative[cbind(
      seq_len(n_penalized),
      n_reactions + seq_len(n_penalized)
    )] <- -1
    blocks <- c(blocks, list(positive, negative))
    lhs <- c(lhs, rep(-Inf, 2L * n_penalized))
    rhs <- c(rhs, rep(0, 2L * n_penalized))
  }

  n_active <- length(active_core)
  core_constraint <- Matrix::Matrix(
    0,
    nrow = n_active,
    ncol = n_reactions + n_penalized,
    sparse = TRUE
  )
  core_constraint[cbind(
    seq_len(n_active),
    match(active_core, colnames(S))
  )] <- 1
  blocks <- c(blocks, list(core_constraint))
  lhs <- c(lhs, rep(scaled_epsilon, n_active))
  rhs <- c(rhs, rep(Inf, n_active))

  auxiliary_upper <- if (n_penalized) {
    pmax(
      abs(scaled_lb[penalized_reactions]),
      abs(scaled_ub[penalized_reactions])
    )
  } else {
    numeric()
  }
  answer <- rc_solve_lp(
    obj = c(rep(0, n_reactions), rep(1, n_penalized)),
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
      status = answer$status, flux = numeric(),
      new_support = character(), objective = NA_real_
    ))
  }
  flux <- answer$solution[seq_len(n_reactions)]
  names(flux) <- colnames(S)
  new_support <- penalized_reactions[
    abs(flux[penalized_reactions]) >= epsilon * (1 - 1e-7)
  ]
  list(
    status = answer$status,
    flux = flux,
    new_support = new_support,
    objective = answer$objective,
    scaling_factor = scaling_factor
  )
}

.rc_fastcore_parent <- function(gem, medium_table = NULL,
                                condition = NULL,
                                forbidden_roles = c(
                                  "demand", "sink",
                                  "artificial_support"
                                ),
                                solver = "highs",
                                time_limit = 300,
                                fastcore_epsilon = 1e-4) {
  parent <- rc_build_full_gem(
    gem,
    medium_table = medium_table,
    condition = condition
  )
  parent <- rc_annotate_reaction_roles(
    parent,
    medium_table = medium_table
  )
  validated <- rc_validate_gem(parent)
  meta <- parent$reaction_meta[
    match(
      validated$reactions,
      as.character(parent$reaction_meta$reaction_id)
    ),
    , drop = FALSE
  ]
  role <- if ("role" %in% colnames(meta)) {
    as.character(meta$role)
  } else {
    rep("unknown", nrow(meta))
  }
  forbidden <- validated$reactions[role %in% forbidden_roles]
  if (length(forbidden)) {
    parent$lb[forbidden] <- 0
    parent$ub[forbidden] <- 0
  }

  validated <- rc_validate_gem(parent)
  feasibility <- rc_solve_lp(
    obj = rep(0, length(validated$reactions)),
    A = validated$S,
    lhs = rep(0, nrow(validated$S)),
    rhs = rep(0, nrow(validated$S)),
    lb = validated$lb,
    ub = validated$ub,
    solver = solver,
    time_limit = time_limit
  )
  if (!identical(feasibility$status, "optimal")) {
    stop(
      "The medium-constrained parent GEM is not feasible: ",
      feasibility$status,
      call. = FALSE
    )
  }

  original_lb <- validated$lb
  original_ub <- validated$ub
  consistent <- .rc_fastcc_consistent_reactions(
    parent,
    solver = solver,
    time_limit = time_limit,
    epsilon = fastcore_epsilon
  )
  inconsistent <- setdiff(validated$reactions, consistent)
  if (length(inconsistent)) {
    parent$lb[inconsistent] <- 0
    parent$ub[inconsistent] <- 0
  }
  parent$fastcore_forbidden_reactions <- forbidden
  parent$fastcc_original_lb <- original_lb
  parent$fastcc_original_ub <- original_ub
  parent$fastcc_consistent_reactions <- consistent
  parent$fastcc_inconsistent_reactions <- inconsistent
  parent
}

.rc_fastcore_complete_direction <- function(parent,
                                             biological_reactions,
                                             selected_support,
                                             targets,
                                             direction,
                                             epsilon,
                                             solver,
                                             time_limit,
                                             max_support_reactions,
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
    validated$S, validated$lb, validated$ub,
    reverse_targets
  )
  remaining <- unique(as.character(targets$reaction_id))
  iteration_rows <- list()
  iteration <- 0L

  local_feasible <- function(reactions, core_ids) {
    local <- .rc_subset_gem(
      list(S = oriented$S, lb = oriented$lb, ub = oriented$ub),
      reactions
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
      lp10 <- .rc_fastcore_lp10(
        oriented$S, oriented$lb, oriented$ub,
        active, penalized, epsilon,
        solver, time_limit,
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
          colnames(oriented$S),
          union(biological_reactions, selected_support)
        )
        one <- .rc_fastcore_lp10(
          oriented$S, oriented$lb, oriented$ub,
          reaction, penalized, epsilon,
          solver, time_limit,
          scaling_factor = scaling_factor
        )
        if (identical(one$status, "optimal")) {
          selected_support <- union(
            selected_support,
            one$new_support
          )
          added <- union(added, one$new_support)
        }
      }
    } else {
      selected_support <- union(selected_support, added)
    }

    if (length(selected_support) > max_support_reactions) {
      stop(
        "FASTCORE support completion exceeded `max_support_reactions`.",
        call. = FALSE
      )
    }

    current_set <- union(biological_reactions, selected_support)
    check <- local_feasible(current_set, remaining)
    remaining <- as.character(
      check$reaction_id[!check$feasible]
    )
    repaired <- setdiff(before, remaining)
    iteration_rows[[iteration]] <- data.frame(
      iteration = iteration,
      target_direction = direction,
      n_targets_before = length(before),
      n_targets_lp7_active = length(active),
      n_targets_repaired = length(repaired),
      n_new_support = length(added),
      lp7_status = lp7$status,
      lp10_status = lp10_status,
      support_objective = objective,
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

.rc_complete_meta_module <- function(gem, reaction_membership,
                                     core_reactions,
                                     sample_id, module_id,
                                     medium_table = NULL,
                                     condition = NULL,
                                     parent_gem = NULL,
                                     target_direction = "both",
                                     solver = "highs",
                                     time_limit = 300,
                                     fastcore_epsilon = 1e-4,
                                     max_support_reactions = 2000,
                                     strict = TRUE) {
  if (!is.finite(fastcore_epsilon) || fastcore_epsilon <= 0) {
    stop("`fastcore_epsilon` must be positive.", call. = FALSE)
  }
  if (!is.finite(max_support_reactions) || max_support_reactions < 0) {
    stop(
      "`max_support_reactions` must be non-negative.",
      call. = FALSE
    )
  }
  required <- c("sample_id", "module_id", "reaction_id")
  if (!is.data.frame(reaction_membership) ||
      !all(required %in% colnames(reaction_membership))) {
    stop(
      paste(
        "`reaction_membership` must contain sample_id,",
        "module_id and reaction_id."
      ),
      call. = FALSE
    )
  }
  if (!is.data.frame(core_reactions) ||
      !all(required %in% colnames(core_reactions))) {
    stop(
      paste(
        "`core_reactions` must contain sample_id,",
        "module_id and reaction_id."
      ),
      call. = FALSE
    )
  }

  parent <- parent_gem %||% .rc_fastcore_parent(
    gem,
    medium_table = medium_table,
    condition = condition,
    solver = solver,
    time_limit = time_limit,
    fastcore_epsilon = fastcore_epsilon
  )
  validated <- rc_validate_gem(parent)
  in_group <-
    as.character(reaction_membership$sample_id) ==
      as.character(sample_id) &
    as.character(reaction_membership$module_id) ==
      as.character(module_id)
  biological <- intersect(
    unique(as.character(
      reaction_membership$reaction_id[in_group]
    )),
    validated$reactions
  )
  if (!length(biological)) {
    stop(
      "No biological reactions found for the requested sample/module.",
      call. = FALSE
    )
  }

  core_group <-
    as.character(core_reactions$sample_id) ==
      as.character(sample_id) &
    as.character(core_reactions$module_id) ==
      as.character(module_id)
  core <- intersect(
    unique(as.character(core_reactions$reaction_id[core_group])),
    biological
  )
  if (!length(core)) {
    stop(
      "No core reactions found for the requested sample/module.",
      call. = FALSE
    )
  }

  direction_model <- parent
  direction_model$lb <- parent$fastcc_original_lb %||% parent$lb
  direction_model$ub <- parent$fastcc_original_ub %||% parent$ub
  target_directions <- rc_prepare_directional_targets(
    direction_model,
    core,
    target_direction = target_direction
  )
  parent_diagnostics <- .rc_directional_feasibility(
    parent,
    target_directions,
    solver = solver,
    time_limit = time_limit,
    flux_threshold = fastcore_epsilon
  )
  parent_feasible_targets <- parent_diagnostics[
    parent_diagnostics$feasible,
    c("reaction_id", "target_direction"),
    drop = FALSE
  ]

  initial <- .rc_subset_gem(parent, biological)
  initial_diagnostics <- .rc_directional_feasibility(
    initial,
    parent_feasible_targets,
    solver = solver,
    time_limit = time_limit,
    flux_threshold = fastcore_epsilon
  )
  names(initial_diagnostics)[
    names(initial_diagnostics) == "feasible"
  ] <- "initial_feasible"
  names(initial_diagnostics)[
    names(initial_diagnostics) == "vmax"
  ] <- "initial_vmax"
  names(initial_diagnostics)[
    names(initial_diagnostics) == "solver_status"
  ] <- "initial_solver_status"
  blocked <- initial_diagnostics[
    !initial_diagnostics$initial_feasible,
    c("reaction_id", "target_direction"),
    drop = FALSE
  ]

  selected_support <- character()
  completion_iterations <- list()
  unresolved <- data.frame(
    reaction_id = character(),
    target_direction = character(),
    stringsAsFactors = FALSE
  )
  for (direction in c("forward", "reverse")) {
    task <- blocked[
      blocked$target_direction == direction,
      , drop = FALSE
    ]
    if (!nrow(task)) next
    completed <- .rc_fastcore_complete_direction(
      parent = parent,
      biological_reactions = biological,
      selected_support = selected_support,
      targets = task,
      direction = direction,
      epsilon = fastcore_epsilon,
      solver = solver,
      time_limit = time_limit,
      max_support_reactions = max_support_reactions
    )
    selected_support <- completed$support
    unresolved <- rbind(unresolved, completed$unresolved)
    if (nrow(completed$iterations)) {
      completion_iterations[[direction]] <- completed$iterations
    }
  }

  final_reactions <- union(biological, selected_support)
  final <- .rc_subset_gem(parent, final_reactions)
  final_diagnostics <- .rc_directional_feasibility(
    final,
    parent_feasible_targets,
    solver = solver,
    time_limit = time_limit,
    flux_threshold = fastcore_epsilon
  )
  names(final_diagnostics)[
    names(final_diagnostics) == "feasible"
  ] <- "final_feasible"
  names(final_diagnostics)[
    names(final_diagnostics) == "vmax"
  ] <- "final_vmax"
  names(final_diagnostics)[
    names(final_diagnostics) == "solver_status"
  ] <- "final_solver_status"

  diagnostics <- merge(
    parent_diagnostics,
    initial_diagnostics,
    by = c("reaction_id", "target_direction"),
    all.x = TRUE,
    sort = FALSE
  )
  diagnostics <- merge(
    diagnostics,
    final_diagnostics,
    by = c("reaction_id", "target_direction"),
    all.x = TRUE,
    sort = FALSE
  )
  diagnostics$completion_status <- ifelse(
    diagnostics$target_direction == "none",
    "no_allowed_direction",
    ifelse(
      !diagnostics$feasible,
      "parent_blocked",
      ifelse(
        diagnostics$initial_feasible %in% TRUE,
        "already_feasible",
        ifelse(
          diagnostics$final_feasible %in% TRUE,
          "fastcore_completed",
          "unresolved"
        )
      )
    )
  )

  failed <- diagnostics$feasible %in% TRUE &
    !(diagnostics$final_feasible %in% TRUE)
  if (isTRUE(strict) && any(failed)) {
    bad <- paste(
      paste(
        diagnostics$reaction_id[failed],
        diagnostics$target_direction[failed],
        sep = ":"
      ),
      collapse = ", "
    )
    stop(
      paste0(
        "FASTCORE meta-module completion failed for ",
        "parent-feasible targets: ", bad
      ),
      call. = FALSE
    )
  }

  meta <- final$reaction_meta
  if (is.null(meta)) {
    meta <- data.frame(
      reaction_id = colnames(final$S),
      stringsAsFactors = FALSE
    )
  }
  meta$biological_meta_module_member <-
    as.character(meta$reaction_id) %in% biological
  meta$fastcore_support <-
    as.character(meta$reaction_id) %in% selected_support
  meta$support_only <- !meta$biological_meta_module_member
  final$reaction_meta <- meta
  final$sample_id <- as.character(sample_id)
  final$grn_module_id <- as.character(module_id)
  final$target_directions <- parent_feasible_targets
  final$closure_diagnostics <- diagnostics
  final$completion_iterations <- if (length(completion_iterations)) {
    do.call(rbind, completion_iterations)
  } else {
    data.frame()
  }

  n_no_direction <- sum(diagnostics$completion_status == "no_allowed_direction")
  n_parent_blocked <- sum(diagnostics$completion_status == "parent_blocked")
  final$target_status <- if (any(failed)) {
    "structurally_infeasible"
  } else if (nrow(diagnostics) > 0L && n_no_direction == nrow(diagnostics)) {
    "no_allowed_direction"
  } else if (n_no_direction > 0L) {
    "partial_no_allowed_direction"
  } else if (nrow(diagnostics) > 0L && n_parent_blocked == nrow(diagnostics)) {
    "parent_blocked"
  } else if (n_parent_blocked > 0L) {
    "partial_parent_blocked"
  } else {
    "ok"
  }
  final$build_params <- list(
    strategy = "meta_module_gem",
    algorithm = "add_only_direction_preserving_fastcore_lp7_lp10",
    n_biological_reactions = length(biological),
    n_fastcore_support_reactions = length(selected_support),
    n_fastcc_consistent_parent_reactions = length(
      parent$fastcc_consistent_reactions %||% colnames(parent$S)
    ),
    n_fastcc_inconsistent_parent_reactions = length(
      parent$fastcc_inconsistent_reactions %||% character()
    ),
    fastcore_epsilon = fastcore_epsilon,
    target_direction = target_direction,
    biological_reactions = biological,
    core_reactions = core,
    forbidden_roles = c(
      "demand", "sink", "artificial_support"
    ),
    strict = strict
  )
  final
}

#' Build one GRN-defined meta-module GEM with add-only FASTCORE completion
.rc_build_meta_module_gem_core <- function(gem, reaction_membership,
                                     core_reactions = NULL,
                                     sample_id, module_id,
                                     medium_table = NULL,
                                     condition = NULL,
                                     target_direction = c(
                                       "both", "forward", "reverse"
                                     ),
                                     solver = "highs",
                                     time_limit = 300,
                                     fastcore_epsilon = 1e-4,
                                     max_support_reactions = 2000,
                                     strict = TRUE) {
  target_direction <- match.arg(target_direction)
  if (is.null(core_reactions)) {
    if (!"is_core" %in% colnames(reaction_membership)) {
      stop(
        paste(
          "Supply `core_reactions` or an `is_core` column",
          "in `reaction_membership`."
        ),
        call. = FALSE
      )
    }
    core_reactions <- reaction_membership[
      reaction_membership$is_core %in% TRUE,
      , drop = FALSE
    ]
  } else if ("is_core" %in% colnames(core_reactions)) {
    core_reactions <- core_reactions[
      core_reactions$is_core %in% TRUE,
      , drop = FALSE
    ]
  }
  .rc_complete_meta_module(
    gem = gem,
    reaction_membership = reaction_membership,
    core_reactions = core_reactions,
    sample_id = sample_id,
    module_id = module_id,
    medium_table = medium_table,
    condition = condition,
    target_direction = target_direction,
    solver = solver,
    time_limit = time_limit,
    fastcore_epsilon = fastcore_epsilon,
    max_support_reactions = max_support_reactions,
    strict = strict
  )
}

rc_build_meta_module_gem <- function(gem, reaction_membership,
                                     core_reactions = NULL, ...) {
  if (!is.null(core_reactions)) {
    core_reactions <- .rc_hard_core_rows(core_reactions)
  }
  .rc_build_meta_module_gem_core(
    gem = gem,
    reaction_membership = reaction_membership,
    core_reactions = core_reactions,
    ...
  )
}
