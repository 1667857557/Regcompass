# Corrected implementation of the resendislab/corda Python CORDA2 flow.

.rc_layer2_corda_options_before_corda2 <- .rc_layer2_corda_options

.rc_layer2_corda_options <- function(model_params = list()) {
  if (!is.list(model_params)) {
    stop("`layer2_args$model_params` must be a list.", call. = FALSE)
  }
  requested <- as.character(model_params$model_completion %||% "fastcore")
  if (length(requested) != 1L || is.na(requested)) {
    stop(
      "`model_completion` must be `fastcore` or `corda2`.",
      call. = FALSE
    )
  }
  is_corda2 <- requested %in% c("corda2", "corda", "corda_like")
  translated <- model_params
  translated$model_completion <- if (is_corda2) "corda" else requested
  translated$corda_gamma <-
    model_params$corda2_penalty_factor %||%
    model_params$corda_penalty_factor %||%
    model_params$corda_gamma %||% 100
  translated$corda_kappa <-
    model_params$corda2_cost_increase %||%
    model_params$corda_cost_increase %||%
    model_params$corda_kappa %||% 1.01
  translated$corda_epsilon <-
    model_params$corda2_target_flux %||%
    model_params$corda_tflux %||%
    model_params$corda_epsilon %||% 1
  translated$corda_n <-
    model_params$corda2_redundancies %||%
    model_params$corda_n %||% 3L
  translated$corda_p <-
    model_params$corda2_support %||%
    model_params$corda_support %||%
    model_params$corda_p %||% 5L
  translated$corda_flux_tolerance <-
    model_params$corda2_flux_tolerance %||%
    model_params$corda_flux_tolerance %||% 1e-8
  answer <- .rc_layer2_corda_options_before_corda2(translated)
  if (!is_corda2) return(answer)
  answer$model_completion <- "corda"
  answer$requested_model_completion <- "corda2"
  answer$penalty_factor <- .rc_corda_scalar(
    translated$corda_gamma,
    "corda2_penalty_factor", 1, Inf
  )
  answer$cost_increase <- .rc_corda_scalar(
    translated$corda_kappa,
    "corda2_cost_increase", 1, Inf
  )
  if (answer$cost_increase <= 1) {
    stop("`corda2_cost_increase` must be greater than 1.", call. = FALSE)
  }
  answer$target_flux <- .rc_corda_scalar(
    translated$corda_epsilon,
    "corda2_target_flux", .Machine$double.eps, Inf
  )
  answer$redundancies <- .rc_corda_integer(
    translated$corda_n,
    "corda2_redundancies", 1L
  )
  answer$support <- .rc_corda_integer(
    translated$corda_p,
    "corda2_support", 1L
  )
  answer$upper_bound <- 1e6
  answer$algorithm <-
    "resendislab_python_CORDA2_corrected_redundant_path_assessment"
  answer$python_reference_commit <-
    "c02e06d50606bf93f23d8f2e6d6ade0e996ca70e"
  answer$python_compatibility <- paste(
    "five directional confidence levels, multiplicative CI redundancy",
    "search, absent-reaction support counting, independent medium testing,",
    "and final free-reaction completion"
  )
  answer$intentional_corrections <- c(
    "maximize remaining medium-confidence directional flux instead of",
    "minimizing a positive target coefficient",
    "block the opposite direction of a reversible target during testing"
  )
  answer
}

.rc_corda2_normalize_split <- function(split, upper_bound = 1e6) {
  answer <- split
  open <- answer$ub > answer$tolerance
  answer$ub[open] <- upper_bound
  answer$algorithm <- "cobra_direction_variables_with_CORDA2_open_bounds"
  answer$corda2_upper_bound <- upper_bound
  answer
}

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
  stats::setNames(value, split$direction_table$variable_id)
}

.rc_corda2_penalties <- function(
    directional_confidence, penalize_medium, penalty_factor) {
  penalty <- stats::setNames(
    rep(0, length(directional_confidence)), names(directional_confidence)
  )
  if (isTRUE(penalize_medium)) {
    penalty[directional_confidence %in% c(1L, 2L)] <- 1
  }
  penalty[directional_confidence == -1L] <- penalty_factor
  penalty
}

.rc_corda2_associated_target <- function(
    engine, target, directional_confidence, options,
    penalize_medium = TRUE, redundancies = TRUE, stage) {
  split <- engine$split
  target <- as.character(target)
  max_iter <- if (isTRUE(redundancies)) options$redundancies else 1L
  penalty <- .rc_corda2_penalties(
    directional_confidence,
    penalize_medium = penalize_medium,
    penalty_factor = options$penalty_factor
  )
  if (!target %in% names(split$ub) || split$ub[[target]] < options$flux_tolerance) {
    return(list(
      result = list(
        target = target,
        status = "target_blocked",
        associated = character(),
        target_flux = 0,
        objective = NA_real_,
        backend = engine$type,
        solver_message = "target upper bound is below solver tolerance",
        opposite_direction_blocked = character(),
        redundancies = 0L,
        n_solves = 0L,
        stage = stage
      ),
      engine = engine,
      impossible = target
    ))
  }
  bounds <- .rc_corda_target_bounds(
    split, target, epsilon = options$target_flux
  )
  needed <- character()
  has_new <- TRUE
  iteration <- 0L
  redundancy_count <- 0L
  status <- "not_run"
  target_flux <- NA_real_
  objective_value <- NA_real_
  backend <- engine$type
  solver_message <- ""
  while (isTRUE(has_new) && iteration < max_iter) {
    iteration <- iteration + 1L
    solved <- .rc_corda_engine_solve(
      engine,
      objective = as.numeric(penalty),
      lower = bounds$lower,
      upper = bounds$upper
    )
    engine <- solved$engine
    answer <- solved$answer
    status <- answer$status
    objective_value <- answer$objective
    backend <- answer$backend
    solver_message <- answer$solver_message %||% ""
    if (!identical(status, "optimal") ||
        length(answer$solution) != ncol(split$S)) {
      break
    }
    flux <- as.numeric(answer$solution)
    names(flux) <- colnames(split$S)
    target_flux <- flux[[target]]
    candidate <- names(directional_confidence)[
      flux[names(directional_confidence)] > options$flux_tolerance &
        directional_confidence %in% c(-1L, 1L, 2L)
    ]
    candidate <- setdiff(candidate, target)
    new <- setdiff(candidate, needed)
    has_new <- length(new) > 0L
    if (isTRUE(redundancies) && has_new) {
      redundancy_count <- redundancy_count + 1L
    }
    weighted_new <- intersect(new, names(penalty)[penalty > 0])
    if (length(weighted_new)) {
      penalty[weighted_new] <-
        penalty[weighted_new] * options$cost_increase
    }
    needed <- union(needed, candidate)
  }
  impossible <- if (identical(status, "optimal")) character() else target
  list(
    result = list(
      target = target,
      status = status,
      associated = sort(unique(needed)),
      target_flux = target_flux,
      objective = objective_value,
      backend = backend,
      solver_message = solver_message,
      opposite_direction_blocked = bounds$opposite_variables,
      redundancies = redundancy_count,
      n_solves = iteration,
      stage = stage
    ),
    engine = engine,
    impossible = impossible
  )
}

.rc_corda2_associated <- function(
    split, targets, directional_confidence, options, solver, time_limit,
    penalize_medium = TRUE, redundancies = TRUE, stage,
    BPPARAM = .rc_corda_task_bpparam()) {
  targets <- intersect(unique(as.character(targets)), colnames(split$S))
  if (!length(targets)) {
    return(list(
      needed = character(),
      results = list(),
      impossible = character(),
      execution = list(
        n_targets = 0L, n_chunks = 0L, workers = 1L,
        solver_runtime = "not_run", n_solves = 0L, n_fallback = 0L
      )
    ))
  }
  workers <- .rc_corda_worker_count(BPPARAM, length(targets))
  n_chunks <- min(workers, length(targets))
  chunk_id <- rep(seq_len(n_chunks), length.out = length(targets))
  chunks <- split(targets, chunk_id)
  run_chunk <- function(chunk) {
    engine <- .rc_corda_new_lp_engine(split, solver, time_limit)
    results <- vector("list", length(chunk))
    impossible <- character()
    for (i in seq_along(chunk)) {
      solved <- .rc_corda2_associated_target(
        engine = engine,
        target = chunk[[i]],
        directional_confidence = directional_confidence,
        options = options,
        penalize_medium = penalize_medium,
        redundancies = redundancies,
        stage = stage
      )
      engine <- solved$engine
      results[[i]] <- solved$result
      impossible <- union(impossible, solved$impossible)
    }
    list(
      results = results,
      impossible = impossible,
      engine = list(
        type = engine$type,
        n_solves = engine$n_solves,
        n_fallback = engine$n_fallback
      )
    )
  }
  parts <- rc_parallel_lapply(
    chunks,
    run_chunk,
    BPPARAM = if (n_chunks > 1L) BPPARAM else FALSE
  )
  results <- unlist(lapply(parts, `[[`, "results"), recursive = FALSE)
  impossible <- unique(unlist(
    lapply(parts, `[[`, "impossible"), use.names = FALSE
  ))
  needed <- unlist(
    lapply(results, `[[`, "associated"), use.names = FALSE
  )
  engine_rows <- do.call(rbind, lapply(parts, function(part) {
    data.frame(
      type = part$engine$type,
      n_solves = part$engine$n_solves,
      n_fallback = part$engine$n_fallback,
      stringsAsFactors = FALSE
    )
  }))
  list(
    needed = as.character(needed),
    results = results,
    impossible = impossible,
    execution = list(
      n_targets = length(targets),
      n_chunks = n_chunks,
      workers = workers,
      task_granularity = "signed_target_with_serial_redundancy_iterations",
      stage_barrier = TRUE,
      persistent_solver = any(engine_rows$type == "highs_persistent_cpp"),
      solver_runtime = paste(unique(engine_rows$type), collapse = ";"),
      n_solves = sum(engine_rows$n_solves),
      n_fallback = sum(engine_rows$n_fallback)
    )
  )
}

.rc_corda2_maximize_targets <- function(
    split, targets, options, solver, time_limit,
    stage = "stage2_independent_medium_flux",
    BPPARAM = .rc_corda_task_bpparam()) {
  targets <- intersect(unique(as.character(targets)), colnames(split$S))
  if (!length(targets)) {
    return(list(feasible = character(), results = list(), execution = list(
      n_targets = 0L, n_chunks = 0L, workers = 1L,
      solver_runtime = "not_run", n_solves = 0L, n_fallback = 0L
    )))
  }
  workers <- .rc_corda_worker_count(BPPARAM, length(targets))
  n_chunks <- min(workers, length(targets))
  chunks <- split(targets, rep(seq_len(n_chunks), length.out = length(targets)))
  run_chunk <- function(chunk) {
    engine <- .rc_corda_new_lp_engine(split, solver, time_limit)
    rows <- vector("list", length(chunk))
    for (i in seq_along(chunk)) {
      target <- chunk[[i]]
      bounds <- .rc_corda_target_bounds(split, target, epsilon = NULL)
      objective <- rep(0, ncol(split$S))
      objective[[bounds$target_index]] <- -1
      solved <- .rc_corda_engine_solve(
        engine, objective, bounds$lower, bounds$upper
      )
      engine <- solved$engine
      answer <- solved$answer
      flux <- if (identical(answer$status, "optimal") &&
                  length(answer$solution) == ncol(split$S)) {
        as.numeric(answer$solution[[bounds$target_index]])
      } else {
        NA_real_
      }
      rows[[i]] <- list(
        target = target,
        status = answer$status,
        associated = character(),
        target_flux = flux,
        objective = answer$objective,
        backend = answer$backend,
        solver_message = answer$solver_message %||% "",
        opposite_direction_blocked = bounds$opposite_variables,
        redundancies = 0L,
        n_solves = 1L,
        stage = stage
      )
    }
    list(
      results = rows,
      engine = list(
        type = engine$type,
        n_solves = engine$n_solves,
        n_fallback = engine$n_fallback
      )
    )
  }
  parts <- rc_parallel_lapply(
    chunks, run_chunk,
    BPPARAM = if (n_chunks > 1L) BPPARAM else FALSE
  )
  results <- unlist(lapply(parts, `[[`, "results"), recursive = FALSE)
  feasible <- vapply(results, function(result) {
    identical(result$status, "optimal") &&
      is.finite(result$target_flux) &&
      result$target_flux > options$target_flux
  }, logical(1))
  engine_rows <- do.call(rbind, lapply(parts, function(part) {
    data.frame(
      type = part$engine$type,
      n_solves = part$engine$n_solves,
      n_fallback = part$engine$n_fallback,
      stringsAsFactors = FALSE
    )
  }))
  list(
    feasible = vapply(results[feasible], `[[`, character(1), "target"),
    results = results,
    execution = list(
      n_targets = length(targets),
      n_chunks = n_chunks,
      workers = workers,
      task_granularity = "signed_target_maximum_flux",
      stage_barrier = TRUE,
      persistent_solver = any(engine_rows$type == "highs_persistent_cpp"),
      solver_runtime = paste(unique(engine_rows$type), collapse = ";"),
      n_solves = sum(engine_rows$n_solves),
      n_fallback = sum(engine_rows$n_fallback)
    )
  )
}

.rc_corda2_results_table <- function(results, split) {
  if (!length(results)) return(.rc_corda_empty_task_table())
  rows <- lapply(results, function(result) {
    target <- as.character(result$target)
    target_row <- split$direction_table[
      split$direction_table$variable_id == target, , drop = FALSE
    ]
    data.frame(
      variable_id = target,
      reaction_id = as.character(target_row$reaction_id[[1L]]),
      direction = as.character(target_row$direction[[1L]]),
      stage = as.character(result$stage),
      replicate = 1L,
      kind = if (grepl("flux$", result$stage)) "feasibility" else "dependency",
      status = as.character(result$status),
      target_flux = as.numeric(result$target_flux),
      objective = as.numeric(result$objective),
      backend = as.character(result$backend),
      solver_message = as.character(result$solver_message %||% ""),
      noise_namespace = "not_applicable_to_corda2",
      opposite_direction_blocked = paste(
        result$opposite_direction_blocked %||% character(), collapse = ";"
      ),
      n_associated = length(result$associated),
      associated = paste(result$associated, collapse = ";"),
      corda2_redundancies = as.integer(result$redundancies %||% 0L),
      corda2_n_solves = as.integer(result$n_solves %||% 0L),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.rc_corda2_reduce_confidence <- function(
    split, directional_confidence, initial_reaction_confidence) {
  reactions <- unique(as.character(split$direction_table$reaction_id))
  reduced <- stats::setNames(rep(-1L, length(reactions)), reactions)
  for (reaction in reactions) {
    variables <- split$direction_table$variable_id[
      split$direction_table$reaction_id == reaction
    ]
    reduced[[reaction]] <- max(directional_confidence[variables])
  }
  label <- as.character(initial_reaction_confidence[reactions])
  label[reduced == 3L] <- "RE"
  stats::setNames(label, reactions)
}

.rc_corda_build_three_stage <- function(
    split, classes, options, solver, time_limit) {
  split <- .rc_corda2_normalize_split(
    split, upper_bound = options$upper_bound
  )
  confidence <- .rc_corda2_directional_confidence(split, classes)
  initial_directional_confidence <- confidence
  inclusion_stage_direction <- stats::setNames(
    rep(NA_character_, length(confidence)), names(confidence)
  )
  inclusion_stage_direction[confidence == 3L] <- "initial_high_confidence"
  execution <- list()
  task_tables <- list()

  stage1_targets <- names(confidence)[confidence == 3L]
  stage1 <- .rc_corda2_associated(
    split, stage1_targets, confidence, options,
    solver = solver, time_limit = time_limit,
    penalize_medium = TRUE, redundancies = TRUE,
    stage = "corda2_stage1_high_associations"
  )
  confidence[stage1$impossible] <- -1L
  stage1_needed <- unique(stage1$needed)
  newly_stage1 <- stage1_needed[confidence[stage1_needed] != 3L]
  confidence[stage1_needed] <- 3L
  inclusion_stage_direction[newly_stage1] <-
    "corda2_stage1_associated_with_high"
  execution$stage1 <- stage1$execution
  task_tables$stage1 <- .rc_corda2_results_table(stage1$results, split)

  stage2_targets <- names(confidence)[confidence %in% c(1L, 2L)]
  stage2 <- .rc_corda2_associated(
    split, stage2_targets, confidence, options,
    solver = solver, time_limit = time_limit,
    penalize_medium = FALSE, redundancies = TRUE,
    stage = "corda2_stage2_medium_absent_support"
  )
  confidence[stage2$impossible] <- -1L
  absent_needed <- stage2$needed[
    confidence[stage2$needed] == -1L
  ]
  absent_count <- table(absent_needed)
  supported_absent <- names(absent_count)[
    as.integer(absent_count) >= options$support
  ]
  newly_absent <- supported_absent[confidence[supported_absent] != 3L]
  confidence[supported_absent] <- 3L
  inclusion_stage_direction[newly_absent] <-
    "corda2_stage2_absent_support_threshold"
  execution$stage2_association <- stage2$execution
  task_tables$stage2_association <-
    .rc_corda2_results_table(stage2$results, split)

  absent_remaining <- names(confidence)[confidence == -1L]
  split_after_absent <- split
  split_after_absent$ub[absent_remaining] <- pmax(
    0, split_after_absent$lb[absent_remaining]
  )
  medium_remaining <- names(confidence)[confidence %in% c(1L, 2L)]
  medium_flux <- .rc_corda2_maximize_targets(
    split_after_absent, medium_remaining, options,
    solver = solver, time_limit = time_limit
  )
  feasible_medium <- unique(medium_flux$feasible)
  newly_medium <- feasible_medium[confidence[feasible_medium] != 3L]
  confidence[feasible_medium] <- 3L
  inclusion_stage_direction[newly_medium] <-
    "corda2_stage2_independent_medium_flux"
  execution$stage2_medium_flux <- medium_flux$execution
  task_tables$stage2_medium_flux <-
    .rc_corda2_results_table(medium_flux$results, split_after_absent)

  split_stage3 <- split_after_absent
  medium_not_included <- names(confidence)[confidence %in% c(1L, 2L)]
  split_stage3$ub[medium_not_included] <- 0
  split_stage3$lb[medium_not_included] <- pmin(
    split_stage3$lb[medium_not_included],
    split_stage3$ub[medium_not_included]
  )
  confidence[confidence == 0L] <- -1L
  stage3_targets <- names(confidence)[confidence == 3L]
  stage3 <- .rc_corda2_associated(
    split_stage3, stage3_targets, confidence, options,
    solver = solver, time_limit = time_limit,
    penalize_medium = FALSE, redundancies = FALSE,
    stage = "corda2_stage3_free_completion"
  )
  confidence[stage3$impossible] <- -1L
  stage3_needed <- unique(stage3$needed)
  newly_stage3 <- stage3_needed[confidence[stage3_needed] != 3L]
  confidence[stage3_needed] <- 3L
  inclusion_stage_direction[newly_stage3] <-
    "corda2_stage3_free_association"
  execution$stage3 <- stage3$execution
  task_tables$stage3 <- .rc_corda2_results_table(stage3$results, split_stage3)

  included_variables <- names(confidence)[confidence == 3L]
  included_reactions <- unique(as.character(
    split$variable_to_reaction[included_variables]
  ))
  initial_reaction_confidence <- classes$initial_confidence
  final_reaction_confidence <- .rc_corda2_reduce_confidence(
    split, confidence, initial_reaction_confidence
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
        initial_directional_confidence[result$associated] == -1L
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
  stage1_reactions <- unique(as.character(
    split$variable_to_reaction[stage1_needed]
  ))
  stage2_absent_reactions <- unique(as.character(
    split$variable_to_reaction[supported_absent]
  ))
  stage2_medium_reactions <- unique(as.character(
    split$variable_to_reaction[feasible_medium]
  ))
  stage3_reactions <- unique(as.character(
    split$variable_to_reaction[stage3_needed]
  ))
  list(
    included = included_reactions,
    included_directional_variables = included_variables,
    final_directional_confidence = confidence,
    initial_directional_confidence = initial_directional_confidence,
    final_confidence = final_reaction_confidence,
    inclusion_stage = inclusion_stage,
    inclusion_stage_direction = inclusion_stage_direction,
    stage1_associated = stage1_reactions,
    stage2_nc_support_pairs = support_pairs,
    stage2_nc_support_count = absent_count,
    stage2_promoted_nc = stage2_absent_reactions,
    stage2_promoted_mc = stage2_medium_reactions,
    stage3_associated_ot = stage3_reactions,
    blocked_after_stage2 = unique(as.character(
      split$variable_to_reaction[absent_remaining]
    )),
    blocked_before_stage3 = unique(as.character(
      split$variable_to_reaction[union(absent_remaining, medium_not_included)]
    )),
    impossible_directional_targets = unique(c(
      stage1$impossible, stage2$impossible, stage3$impossible
    )),
    redundancies = stats::setNames(
      vapply(c(stage1$results, stage2$results), function(result) {
        as.integer(result$redundancies %||% 0L)
      }, integer(1)),
      vapply(c(stage1$results, stage2$results), `[[`, character(1), "target")
    ),
    task_diagnostics = .rc_bind_frames_fill(task_tables),
    execution = execution,
    algorithm =
      "resendislab_python_CORDA2_corrected_redundant_path_assessment",
    python_reference_commit =
      "c02e06d50606bf93f23d8f2e6d6ade0e996ca70e",
    stage_update_policy = "python_build_stage_barriers",
    corrected_python_defects = c(
      "remaining medium confidence flux is maximized",
      "opposite reversible target direction is blocked"
    )
  )
}

.rc_validate_corda_union_model_before_corda2 <-
  .rc_validate_corda_union_model

.rc_validate_corda_union_model <- function(model, cell_type) {
  algorithm <- as.character(model$build_params$algorithm %||% "")
  if (!identical(
    algorithm,
    "resendislab_python_CORDA2_corrected_redundant_path_assessment"
  )) {
    return(.rc_validate_corda_union_model_before_corda2(model, cell_type))
  }
  copy <- model
  copy$build_params$algorithm <-
    "Schultz_Qutub_CORDA_2016_three_stage_dependency_assessment"
  copy$build_params$stage_update_policy <-
    "barrier_then_union_order_independent"
  .rc_validate_corda_union_model_before_corda2(copy, cell_type)
  if (!identical(
    as.character(model$build_params$python_reference_commit),
    "c02e06d50606bf93f23d8f2e6d6ade0e996ca70e"
  )) {
    stop("CORDA2 Python reference commit is missing.", call. = FALSE)
  }
  invisible(TRUE)
}

.rc_complete_celltype_medium_corda_gem_before_corda2 <-
  .rc_complete_celltype_medium_corda_gem

.rc_complete_celltype_medium_corda_gem <- function(...) {
  model <- .rc_complete_celltype_medium_corda_gem_before_corda2(...)
  reconstruction <- model$corda_reconstruction
  model$build_params$strategy <- "celltype_medium_corrected_python_corda2"
  model$build_params$algorithm <- reconstruction$algorithm
  model$build_params$completion_stage <-
    "corrected_python_CORDA2_after_confidence_mapping"
  model$build_params$stage_update_policy <- reconstruction$stage_update_policy
  model$build_params$python_reference_commit <-
    reconstruction$python_reference_commit
  model$build_params$corrected_python_defects <-
    reconstruction$corrected_python_defects
  model$build_params$corda2_redundancies <-
    model$build_params$corda_options$redundancies %||%
    model$build_params$corda_options$n
  model$build_params$corda2_support <-
    model$build_params$corda_options$support %||%
    model$build_params$corda_options$p
  model$build_params$corda2_penalty_factor <-
    model$build_params$corda_options$penalty_factor %||%
    model$build_params$corda_options$gamma
  model$build_params$corda2_cost_increase <-
    model$build_params$corda_options$cost_increase %||%
    model$build_params$corda_options$kappa
  model$build_params$corda2_target_flux <-
    model$build_params$corda_options$target_flux %||%
    model$build_params$corda_options$epsilon
  model$build_params$corda2_upper_bound <-
    model$build_params$corda_options$upper_bound %||% 1e6
  model$corda2_contract <- list(
    implementation = "corrected resendislab/corda Python CORDA2",
    reference_repository = "resendislab/corda",
    reference_commit = reconstruction$python_reference_commit,
    confidence_levels = c(absent = -1L, unknown = 0L,
                          low = 1L, medium = 2L, high = 3L),
    redundant_path_cost_increase =
      model$build_params$corda2_cost_increase,
    absent_penalty_factor =
      model$build_params$corda2_penalty_factor,
    support_threshold = model$build_params$corda2_support,
    maximum_redundant_paths = model$build_params$corda2_redundancies,
    target_flux = model$build_params$corda2_target_flux,
    corrections = reconstruction$corrected_python_defects
  )
  model
}

.rc_complete_celltype_medium_corda_like_gem <-
  .rc_complete_celltype_medium_corda_gem
