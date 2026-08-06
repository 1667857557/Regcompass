# Directional bounds used by the original MATLAB CORDA2 decomposition.

.rc_corda_target_bounds <- function(split, target, epsilon = NULL) {
  target <- as.character(target)
  if (length(target) != 1L || is.na(target) ||
      !target %in% split$direction_table$variable_id) {
    stop("CORDA2 target variable is not present in the split model.",
         call. = FALSE)
  }
  lower <- split$lb
  upper <- split$ub
  closed <- .rc_corda2_close_opposite(split, target, lower, upper)
  lower <- closed$lower
  upper <- closed$upper
  if (!is.null(epsilon)) {
    epsilon <- as.numeric(epsilon)
    if (length(epsilon) != 1L || !is.finite(epsilon) || epsilon < 0) {
      stop("CORDA2 target flux must be one non-negative finite number.",
           call. = FALSE)
    }
    lower[[target]] <- epsilon
    upper[[target]] <- epsilon
  }
  row <- split$direction_table[
    split$direction_table$variable_id == target, , drop = FALSE
  ]
  list(
    lower = lower,
    upper = upper,
    target = target,
    target_reaction = as.character(row$reaction_id[[1L]]),
    opposite_variables = closed$opposite,
    opposite_direction_blocked = closed$opposite,
    target_index = match(target, names(upper)),
    source_semantics = paste(
      "original CORDA2 closes the opposite directional copy and fixes the",
      "tested target at its val-or-percentage constraint"
    )
  )
}

