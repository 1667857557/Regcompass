# Exact Python CORDA2 build-state transitions.

.rc_corda2_minimize_medium_targets <- function(
    engine, split, targets, directional_confidence, options,
    stage = "corda2_stage2_independent_medium_minimization") {
  targets <- as.character(targets)
  if (anyNA(targets) || any(!targets %in% colnames(split$S))) {
    stop("CORDA2 medium variable is missing from the solver model.",
         call. = FALSE)
  }
  results <- vector("list", length(targets))
  promoted <- character()

  for (i in seq_along(targets)) {
    target <- targets[[i]]
    objective <- stats::setNames(
      rep(0, ncol(split$S)), colnames(split$S)
    )
    objective[[target]] <- 1
    solved <- .rc_corda_engine_solve(
      engine,
      objective = as.numeric(objective),
      lower = split$lb,
      upper = split$ub
    )
    engine <- solved$engine
    answer <- solved$answer
    target_flux <- if (identical(answer$status, "optimal") &&
                       length(answer$solution) == ncol(split$S)) {
      as.numeric(answer$solution[[match(target, colnames(split$S))]])
    } else {
      NA_real_
    }
    if (identical(answer$status, "optimal") &&
        is.finite(answer$objective) &&
        answer$objective > options$target_flux) {
      directional_confidence[[target]] <- 3L
      promoted <- c(promoted, target)
    }
    results[[i]] <- .rc_corda2_target_result(
      split, target, stage, "positive_coefficient_minimization",
      answer$status,
      target_flux = target_flux,
      objective = answer$objective,
      backend = answer$backend,
      solver_message = answer$solver_message %||% "",
      n_solves = 1L
    )
  }

  list(
    engine = engine,
    confidence = directional_confidence,
    promoted = promoted,
    results = results,
    execution = list(
      n_targets = length(targets),
      n_chunks = if (length(targets)) 1L else 0L,
      workers = 1L,
      task_granularity = "python_serial_variable_order_positive_minimization",
      stage_barrier = FALSE,
      target_parallelism = FALSE,
      persistent_solver = identical(engine$type, "highs_persistent_cpp"),
      solver_runtime = engine$type,
      n_solves = engine$n_solves,
      n_fallback = engine$n_fallback
    )
  )
}

.rc_corda2_results_table <- function(results, split = NULL) {
  if (!length(results)) return(.rc_corda_empty_task_table())
  .rc_corda_results_table(results)
}

# Exact equivalent of Python CORDA.__reduce_conf(): one numeric confidence per
# reaction, equal to max(forward confidence, reverse confidence).
.rc_corda2_reduce_confidence <- function(split, directional_confidence) {
  reactions <- split$reaction_order
  reduced <- stats::setNames(rep(-1L, length(reactions)), reactions)
  for (reaction in reactions) {
    variables <- split$direction_table$variable_id[
      split$direction_table$reaction_id == reaction
    ]
    reduced[[reaction]] <- max(directional_confidence[variables])
  }
  reduced
}

.rc_corda2_count_in_order <- function(value) {
  value <- as.character(value)
  if (!length(value)) return(stats::setNames(integer(), character()))
  keys <- unique(value)
  stats::setNames(vapply(keys, function(key) sum(value == key), integer(1)), keys)
}

.rc_corda_build_three_stage <- function(
    split, classes, options, solver, time_limit) {
  confidence <- .rc_corda2_directional_confidence(split, classes)
  initial_directional_confidence <- confidence
  inclusion_stage_direction <- stats::setNames(
    rep(NA_character_, length(confidence)), names(confidence)
  )
  inclusion_stage_direction[confidence == 3L] <- "initial_high_confidence"
  engine <- .rc_corda_new_lp_engine(split, solver, time_limit)
  execution <- list()
  task_tables <- list()
  all_impossible <- character()
  redundancy_values <- stats::setNames(
    integer(length(confidence)), names(confidence)
  )

  # Python build iteration 1.
  stage1_targets <- names(confidence)[confidence == 3L]
  stage1 <- .rc_corda2_associated(
    engine, split, stage1_targets, confidence, options,
    penalize_medium = TRUE, redundancies = TRUE,
    stage = "corda2_stage1_high_associations"
  )
  engine <- stage1$engine
  confidence <- stage1$confidence
  all_impossible <- c(all_impossible, stage1$impossible)
  redundancy_values[names(stage1$redundancies)] <- stage1$redundancies
  stage1_needed <- unique(stage1$needed)
  newly_stage1 <- stage1_needed[confidence[stage1_needed] != 3L]
  confidence[stage1_needed] <- 3L
  confidence_after_stage1 <- confidence
  inclusion_stage_direction[newly_stage1] <-
    "corda2_stage1_associated_with_high"
  execution$stage1 <- stage1$execution
  task_tables$stage1 <- .rc_corda2_results_table(stage1$results)

  # Python build iteration 2: absent support from low/medium targets.
  stage2_targets <- names(confidence)[confidence %in% c(1L, 2L)]
  stage2 <- .rc_corda2_associated(
    engine, split, stage2_targets, confidence, options,
    penalize_medium = FALSE, redundancies = TRUE,
    stage = "corda2_stage2_medium_absent_support"
  )
  engine <- stage2$engine
  confidence <- stage2$confidence
  all_impossible <- c(all_impossible, stage2$impossible)
  redundancy_values[names(stage2$redundancies)] <- stage2$redundancies
  confidence_after_stage2_association <- confidence
  absent_needed <- stage2$needed[
    confidence_after_stage2_association[stage2$needed] == -1L
  ]
  absent_count <- .rc_corda2_count_in_order(absent_needed)
  supported_absent <- names(absent_count)[
    as.integer(absent_count) >= options$support
  ]
  newly_absent <- supported_absent[confidence[supported_absent] != 3L]
  confidence[supported_absent] <- 3L
  confidence_after_stage2_support <- confidence
  inclusion_stage_direction[newly_absent] <-
    "corda2_stage2_absent_support_threshold"
  execution$stage2_association <- stage2$execution
  task_tables$stage2_association <- .rc_corda2_results_table(stage2$results)

  # Exact Python blocking of all remaining confidence -1 variables.
  split_after_absent <- split
  absent_remaining <- names(confidence)[confidence == -1L]
  split_after_absent$ub[absent_remaining] <- pmax(
    0, split_after_absent$lb[absent_remaining]
  )

  # Exact source behavior: minimize +1 * v and compare objective > tflux.
  medium_targets <- names(confidence)[confidence %in% c(1L, 2L)]
  stage2_medium <- .rc_corda2_minimize_medium_targets(
    engine, split_after_absent, medium_targets, confidence, options
  )
  engine <- stage2_medium$engine
  confidence <- stage2_medium$confidence
  confidence_after_stage2_medium <- confidence
  feasible_medium <- stage2_medium$promoted
  inclusion_stage_direction[feasible_medium] <-
    "corda2_stage2_positive_coefficient_minimization"
  execution$stage2_medium <- stage2_medium$execution
  task_tables$stage2_medium <-
    .rc_corda2_results_table(stage2_medium$results)

  # Python build iteration 3.
  split_stage3 <- split_after_absent
  for (variable in names(confidence)) {
    if (confidence[[variable]] %in% c(1L, 2L)) {
      if (split_stage3$lb[[variable]] > 0) {
        stop(
          "Exact Python CORDA2 cannot set the stage-3 upper bound to zero ",
          "when the variable has a positive lower bound: ", variable,
          call. = FALSE
        )
      }
      split_stage3$ub[[variable]] <- 0
    } else if (identical(confidence[[variable]], 0L)) {
      confidence[[variable]] <- -1L
    }
  }
  stage3_targets <- names(confidence)[confidence == 3L]
  stage3 <- .rc_corda2_associated(
    engine, split_stage3, stage3_targets, confidence, options,
    penalize_medium = FALSE, redundancies = FALSE,
    stage = "corda2_stage3_free_completion"
  )
  engine <- stage3$engine
  confidence <- stage3$confidence
  all_impossible <- c(all_impossible, stage3$impossible)
  stage3_needed <- unique(stage3$needed)
  newly_stage3 <- stage3_needed[confidence[stage3_needed] != 3L]
  confidence[stage3_needed] <- 3L
  inclusion_stage_direction[newly_stage3] <-
    "corda2_stage3_free_association"
  execution$stage3 <- stage3$execution
  task_tables$stage3 <- .rc_corda2_results_table(stage3$results)

  impossible <- sort(unique(all_impossible), method = "radix")
  redundancy_values <- redundancy_values[
    confidence[names(redundancy_values)] == 3L &
      !names(redundancy_values) %in% impossible
  ]
  included_variables <- names(confidence)[confidence == 3L]
  included_reactions <- unique(as.character(
    split$variable_to_reaction[included_variables]
  ))
  initial_reaction_confidence <- .rc_corda2_reduce_confidence(
    split, initial_directional_confidence
  )
  final_reaction_confidence <- .rc_corda2_reduce_confidence(
    split, confidence
  )
  final_reaction_status <- stats::setNames(
    ifelse(final_reaction_confidence == 3L, "included", "excluded"),
    names(final_reaction_confidence)
  )
  inclusion_stage <- stats::setNames(
    rep(NA_character_, length(initial_reaction_confidence)),
    names(initial_reaction_confidence)
  )
  for (reaction in names(inclusion_stage)) {
    variables <- split$direction_table$variable_id[
      split$direction_table$reaction_id == reaction &
        confidence[split$direction_table$variable_id] == 3L
    ]
    stages <- inclusion_stage_direction[variables]
    stages <- stages[!is.na(stages) & nzchar(stages)]
    if (length(stages)) inclusion_stage[[reaction]] <- stages[[1L]]
  }

  support_pairs <- data.frame(
    target_variable = character(),
    absent_variable = character(),
    stringsAsFactors = FALSE
  )
  if (length(stage2$results)) {
    pair_rows <- lapply(stage2$results, function(result) {
      absent <- result$associated[
        confidence_after_stage2_association[result$associated] == -1L
      ]
      if (!length(absent)) return(NULL)
      data.frame(
        target_variable = rep(result$target, length(absent)),
        absent_variable = absent,
        stringsAsFactors = FALSE
      )
    })
    pair_rows <- pair_rows[!vapply(pair_rows, is.null, logical(1))]
    if (length(pair_rows)) support_pairs <- do.call(rbind, pair_rows)
  }

  list(
    included = included_reactions,
    included_directional_variables = included_variables,
    initial_reaction_confidence = initial_reaction_confidence,
    final_reaction_confidence = final_reaction_confidence,
    final_reaction_status = final_reaction_status,
    final_confidence = final_reaction_confidence,
    final_directional_confidence = confidence,
    initial_directional_confidence = initial_directional_confidence,
    confidence_after_stage1 = confidence_after_stage1,
    confidence_after_stage2_association =
      confidence_after_stage2_association,
    confidence_after_stage2_support = confidence_after_stage2_support,
    confidence_after_stage2_medium = confidence_after_stage2_medium,
    inclusion_stage = inclusion_stage,
    inclusion_stage_direction = inclusion_stage_direction,
    stage1_associated = unique(as.character(
      split$variable_to_reaction[stage1_needed]
    )),
    stage2_nc_support_pairs = support_pairs,
    stage2_nc_support_count = absent_count,
    stage2_promoted_nc = unique(as.character(
      split$variable_to_reaction[supported_absent]
    )),
    stage2_promoted_mc = unique(as.character(
      split$variable_to_reaction[feasible_medium]
    )),
    stage3_associated_ot = unique(as.character(
      split$variable_to_reaction[stage3_needed]
    )),
    blocked_after_stage2 = unique(as.character(
      split$variable_to_reaction[absent_remaining]
    )),
    blocked_before_stage3 = unique(as.character(
      split$variable_to_reaction[c(
        absent_remaining,
        names(confidence)[confidence %in% c(1L, 2L)]
      )]
    )),
    impossible_directional_targets = impossible,
    redundancies = redundancy_values,
    task_diagnostics = .rc_bind_frames_fill(task_tables),
    execution = execution,
    algorithm = "resendislab_python_CORDA2_c02e06d_exact_semantics",
    python_reference_commit =
      "c02e06d50606bf93f23d8f2e6d6ade0e996ca70e",
    stage_update_policy = "python_serial_mutation_order",
    source_semantics = c(
      "fixed UPPER=1e6, CI=1.01 and tflux=1",
      "both forward and reverse variables exist for every reaction",
      "target assessment changes only the target variable bounds",
      "penalties for both directions use the forward variable confidence",
      "targets are processed serially with one persistent solver state",
      "remaining medium variables use positive-coefficient minimization"
    )
  )
}