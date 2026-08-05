# Helpers matching the original MATLAB CORDA2.m model decomposition and bounds.

.rc_corda2_split_original <- function(gem) {
  validated <- rc_validate_gem(gem)
  reactions <- validated$reactions
  reaction_lb <- as.numeric(validated$lb)
  reaction_ub <- as.numeric(validated$ub)
  names(reaction_lb) <- names(reaction_ub) <- reactions
  reverse_index <- which(reaction_lb < 0 & reaction_ub >= 0)
  reverse_only <- which(reaction_ub < 0)
  if (length(reverse_only)) {
    stop(
      "Original CORDA2.m does not support reactions with an upper bound below zero: ",
      paste(utils::head(reactions[reverse_only], 10L), collapse = ", "),
      call. = FALSE
    )
  }

  forward_ids <- reactions
  reverse_ids <- paste0(reactions[reverse_index], "_CORDA_rev_rxn")
  variable_ids <- c(forward_ids, reverse_ids)
  S <- cbind(
    validated$S,
    if (length(reverse_index)) -validated$S[, reverse_index, drop = FALSE]
    else Matrix::Matrix(0, nrow(validated$S), 0, sparse = TRUE)
  )
  colnames(S) <- variable_ids
  lower <- stats::setNames(rep(0, length(variable_ids)), variable_ids)
  upper <- stats::setNames(c(reaction_ub, -reaction_lb[reverse_index]), variable_ids)
  direction_table <- rbind(
    data.frame(
      variable_id = forward_ids,
      reaction_id = reactions,
      direction = "forward",
      original_index = seq_along(reactions),
      stringsAsFactors = FALSE
    ),
    data.frame(
      variable_id = reverse_ids,
      reaction_id = reactions[reverse_index],
      direction = "reverse",
      original_index = reverse_index,
      stringsAsFactors = FALSE
    )
  )
  rownames(direction_table) <- NULL
  list(
    S = .rc_as_dgCMatrix(S),
    lb = lower,
    ub = upper,
    direction_table = direction_table,
    reaction_order = reactions,
    variable_order = variable_ids,
    variable_to_reaction = stats::setNames(
      direction_table$reaction_id, direction_table$variable_id
    ),
    variable_to_direction = stats::setNames(
      direction_table$direction, direction_table$variable_id
    ),
    original_reaction_lb = reaction_lb,
    original_reaction_ub = reaction_ub,
    reversible_reactions = reactions[reverse_index],
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

.rc_corda2_maximize_target <- function(
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
    output$lb[[reaction]] <- if (keep_reverse) max(parent_lb, -1000) else 0
    output$ub[[reaction]] <- if (keep_forward) min(parent_ub, 1000) else 0
  }
  output
}
