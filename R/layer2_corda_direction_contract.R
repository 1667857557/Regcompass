# Exact target-bound semantics of resendislab/corda CORDA.associated().

.rc_corda_target_bounds <- function(split, target, epsilon = NULL) {
  target <- as.character(target)
  if (length(target) != 1L || is.na(target) ||
      !target %in% split$direction_table$variable_id) {
    stop("CORDA2 target variable is not present in the split model.",
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
  target_index <- match(target, names(upper))
  if (!is.null(epsilon)) {
    epsilon <- as.numeric(epsilon)
    if (length(epsilon) != 1L || !is.finite(epsilon) || epsilon <= 0) {
      stop("CORDA2 target flux must be one positive finite number.",
           call. = FALSE)
    }
    new_lower <- max(lower[[target_index]], epsilon)
    # Python assigns `va.lb` before `va.ub = UPPER`. Optlang rejects the
    # transient state when the new lower bound exceeds the current upper bound.
    if (new_lower > upper[[target_index]]) {
      stop(
        "Exact Python CORDA2 target-bound assignment would set lower bound ",
        new_lower, " above the current upper bound ",
        upper[[target_index]], " for `", target, "`.",
        call. = FALSE
      )
    }
    lower[[target_index]] <- new_lower
    upper[[target_index]] <- split$upper_bound
  }
  list(
    lower = lower,
    upper = upper,
    target = target,
    target_reaction = reaction,
    opposite_variables = opposite,
    opposite_direction_blocked = character(),
    target_index = target_index,
    source_semantics = paste(
      "only the target variable bounds are changed; the opposite reversible",
      "variable remains available exactly as in Python CORDA.associated"
    )
  )
}

.rc_corda_results_table <- function(results) {
  if (!length(results)) return(data.frame())
  rows <- lapply(results, function(x) {
    data.frame(
      variable_id = as.character(x$target),
      reaction_id = as.character(x$reaction_id),
      direction = as.character(x$direction),
      stage = as.character(x$stage),
      replicate = 1L,
      kind = as.character(x$kind),
      status = as.character(x$status),
      target_flux = as.numeric(x$target_flux),
      objective = as.numeric(x$objective),
      backend = as.character(x$backend),
      solver_message = as.character(x$solver_message %||% ""),
      noise_namespace = "not_applicable_to_python_corda2",
      opposite_direction_blocked = "",
      n_associated = length(x$associated),
      associated = paste(x$associated, collapse = ";"),
      corda2_redundancies = as.integer(x$redundancies %||% 0L),
      corda2_n_solves = as.integer(x$n_solves %||% 0L),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
