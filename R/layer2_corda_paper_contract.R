# Helpers matching the original MATLAB CORDA2.m model decomposition and bounds.

.rc_corda2_split_original <- function(gem) {
  validated <- rc_validate_gem(gem)
  reactions <- as.character(validated$reactions)
  reaction_lb <- validated$lb[reactions]
  reaction_ub <- validated$ub[reactions]
  if (length(reaction_lb) != length(reactions) ||
      length(reaction_ub) != length(reactions)) {
    stop(
      "CORDA2 input bounds do not align with the reaction order: reactions=",
      length(reactions), ", lb=", length(reaction_lb),
      ", ub=", length(reaction_ub), ".",
      call. = FALSE
    )
  }
  reaction_lb <- stats::setNames(as.numeric(reaction_lb), reactions)
  reaction_ub <- stats::setNames(as.numeric(reaction_ub), reactions)
  reverse_reactions <- reactions[reaction_lb < 0 & reaction_ub >= 0]
  reverse_only <- reactions[reaction_ub < 0]
  if (length(reverse_only)) {
    stop(
      "Original CORDA2.m does not support reactions with an upper bound below zero: ",
      paste(utils::head(reverse_only, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  forward_ids <- reactions
  reverse_ids <- if (length(reverse_reactions)) {
    paste0(reverse_reactions, "_CORDA_rev_rxn")
  } else {
    character()
  }
  variable_ids <- c(forward_ids, reverse_ids)
  triplet <- Matrix::summary(.rc_as_dgCMatrix(validated$S[, reactions, drop = FALSE]))
  reverse_index <- match(reverse_reactions, reactions)
  reverse_entry <- triplet$j %in% reverse_index
  reverse_position <- match(triplet$j[reverse_entry], reverse_index)
  S <- Matrix::sparseMatrix(
    i = c(triplet$i, triplet$i[reverse_entry]),
    j = c(triplet$j, length(reactions) + reverse_position),
    x = c(triplet$x, -triplet$x[reverse_entry]),
    dims = c(nrow(validated$S), length(variable_ids)),
    dimnames = list(rownames(validated$S), variable_ids),
    giveCsparse = TRUE
  )
  upper_values <- c(
    unname(reaction_ub),
    -unname(reaction_lb[reverse_reactions])
  )
  if (length(upper_values) != length(variable_ids)) {
    stop(
      "CORDA2 directional bounds do not align with variable identifiers: variables=",
      length(variable_ids), ", bounds=", length(upper_values),
      ", reversible=", length(reverse_reactions), ".",
      call. = FALSE
    )
  }
  lower <- stats::setNames(rep(0, length(variable_ids)), variable_ids)
  upper <- stats::setNames(upper_values, variable_ids)
  forward_table <- data.frame(
    variable_id = forward_ids,
    reaction_id = reactions,
    direction = rep("forward", length(forward_ids)),
    original_index = seq_along(reactions),
    stringsAsFactors = FALSE
  )
  reverse_table <- if (length(reverse_ids)) {
    data.frame(
      variable_id = reverse_ids,
      reaction_id = reverse_reactions,
      direction = rep("reverse", length(reverse_ids)),
      original_index = reverse_index,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      variable_id = character(),
      reaction_id = character(),
      direction = character(),
      original_index = integer(),
      stringsAsFactors = FALSE
    )
  }
  direction_table <- rbind(forward_table, reverse_table)
  rownames(direction_table) <- NULL
  if (nrow(direction_table) != length(variable_ids)) {
    stop(
      "CORDA2 direction table does not align with variable identifiers: rows=",
      nrow(direction_table), ", variables=", length(variable_ids), ".",
      call. = FALSE
    )
  }
  list(
    S = .rc_as_dgCMatrix(S),
    lb = lower,
    ub = upper,
    direction_table = direction_table,
    reaction_order = reactions,
    variable_order = variable_ids,
    variable_to_reaction = stats::setNames(
      as.character(direction_table$reaction_id), variable_ids
    ),
    variable_to_direction = stats::setNames(
      as.character(direction_table$direction), variable_ids
    ),
    original_reaction_lb = reaction_lb,
    original_reaction_ub = reaction_ub,
    reversible_reactions = reverse_reactions,
    tolerance = 1e-7,
    upper_bound = max(c(upper, 1), na.rm = TRUE),
    algorithm = "original_CORDA2_irreversible_decomposition"
  )
}

.rc_corda2_opposite_variable <- function(split, target) {
  row <- split$direction_table[
    split$direction_table$variable_id == target, , drop = FALSE
  ]
  if (nrow(row) != 1L) {
    stop("CORDA2 target variable is not present in the split model.",
         call. = FALSE)
  }
  split$direction_table$variable_id[
    split$direction_table$reaction_id == row$reaction_id[[1L]] &
      split$direction_table$variable_id != target
  ]
}

.rc_corda2_close_opposite <- function(split, target, lower, upper) {
  opposite <- .rc_corda2_opposite_variable(split, target)
  if (length(opposite)) {
    lower[opposite] <- 0
    upper[opposite] <- 0
  }
  list(lower = lower, upper = upper, opposite = opposite)
}

.rc_corda2_maximize_target_core <- function(
    engine, split, target, lower = split$lb, upper = split$ub) {
  closed <- .rc_corda2_close_opposite(split, target, lower, upper)
  objective <- stats::setNames(rep(0, ncol(split$S)), colnames(split$S))
  objective[[target]] <- -1
  solved <- .rc_corda_engine_solve(
    engine,
    objective = as.numeric(objective),
    lower = closed$lower,
    upper = closed$upper
  )
  answer <- solved$answer
  flux <- if (identical(answer$status, "optimal") &&
              length(answer$solution) == ncol(split$S)) {
    as.numeric(answer$solution[[match(target, colnames(split$S))]])
  } else {
    NA_real_
  }
  list(
    engine = solved$engine,
    answer = answer,
    vmax = flux,
    lower = closed$lower,
    upper = closed$upper,
    opposite = closed$opposite
  )
}

.rc_corda2_constrain_target <- function(
    engine, split, target, options,
    lower = split$lb, upper = split$ub) {
  maximum <- .rc_corda2_maximize_target(
    engine, split, target, lower = lower, upper = upper
  )
  required <- NA_real_
  if (identical(maximum$answer$status, "optimal") &&
      is.finite(maximum$vmax)) {
    required <- if (identical(options$constrainby, "val")) {
      min(options$constraint, maximum$vmax)
    } else {
      0.01 * options$constraint * maximum$vmax
    }
    maximum$lower[[target]] <- required
    maximum$upper[[target]] <- required
  }
  maximum$required_flux <- required
  maximum
}

.rc_corda2_apply_direction_bounds <- function(parent, included_variables, split) {
  validated <- rc_validate_gem(parent)
  included_variables <- unique(as.character(included_variables))
  selected_reactions <- unique(as.character(
    split$variable_to_reaction[included_variables]
  ))
  selected_reactions <- validated$reactions[
    validated$reactions %in% selected_reactions
  ]
  output <- .rc_subset_gem(parent, selected_reactions)
  for (reaction in colnames(output$S)) {
    variables <- split$direction_table[
      split$direction_table$reaction_id == reaction, , drop = FALSE
    ]
    forward <- variables$variable_id[variables$direction == "forward"]
    reverse <- variables$variable_id[variables$direction == "reverse"]
    keep_forward <- length(forward) && forward %in% included_variables
    keep_reverse <- length(reverse) && reverse %in% included_variables
    parent_lb <- validated$lb[[reaction]]
    parent_ub <- validated$ub[[reaction]]
    output$lb[[reaction]] <- if (keep_reverse) {
      max(parent_lb, -1000)
    } else if (keep_forward && parent_lb > 0) {
      parent_lb
    } else {
      0
    }
    output$ub[[reaction]] <- if (keep_forward) min(parent_ub, 1000) else 0
  }
  output
}

# Progress-aware entry point; the algorithm remains in the core above.
.rc_corda2_maximize_target <- function(...) {
  task <- .rc_layer2_progress_state$current_task
  if (!is.null(task) && identical(task$route, "corda2") &&
      !isTRUE(.rc_layer2_progress_state$inside_dependency)) {
    .rc_layer2_algorithm_once(
      "corda2_step2_2", "corda2_step2_2_MC_feasibility", 6L,
      "promoting frequent NC dependencies and testing MC feasibility"
    )
  }
  do.call(.rc_corda2_maximize_target_core, list(...))
}
