# Stage-barrier parallel execution for original MATLAB CORDA2 semantics.
#
# Scheduling changes only. Every directional target receives the same split
# model, confidence snapshot, objective, bounds and stopping rules as the serial
# implementation. State mutations remain at the original CORDA2 stage barriers.

.rc_corda_build_three_stage_core_serial <- .rc_corda_build_three_stage_core

.rc_corda_stage_backend <- function(BPPARAM) {
  if (identical(BPPARAM, FALSE) || is.null(BPPARAM) ||
      !requireNamespace("BiocParallel", quietly = TRUE) ||
      !methods::is(BPPARAM, "BiocParallelParam")) return("auto")
  if (methods::is(BPPARAM, "SnowParam")) return("snow")
  if (methods::is(BPPARAM, "MulticoreParam")) return("multicore")
  "auto"
}

.rc_corda_stage_param <- function(n_targets) {
  n_targets <- max(0L, as.integer(n_targets[[1L]]))
  if (n_targets <= 1L ||
      !isTRUE(.rc_layer2_parallel_context$active) ||
      !isTRUE(.rc_layer2_parallel_context$parallel)) return(FALSE)
  template <- .rc_layer2_parallel_context$BPPARAM
  requested <- if (identical(template, FALSE)) {
    1L
  } else if (is.null(template)) {
    max(1L, as.integer(rc_available_workers(default = 1L)))
  } else {
    .rc_corda_pool_workers(template)
  }
  workers <- min(n_targets, max(1L, as.integer(requested)))
  if (workers <= 1L) return(FALSE)
  param <- rc_default_bpparam(
    workers = workers,
    backend = .rc_corda_stage_backend(template)
  )
  if (is.null(param)) return(FALSE)
  param <- .rc_layer2_tune_task_bpparam(param, n_targets)
  attr(param, "regcompass_corda2_stage_workers") <- workers
  attr(param, "regcompass_corda2_stage_targets") <- n_targets
  param
}

.rc_corda_stage_label <- function(stage) {
  switch(
    as.character(stage),
    corda2_step1_HC_dependencies = "CORDA2 Step 1 HC dependencies",
    corda2_step2_1_MC_NC_dependencies = "CORDA2 Step 2.1 MC-to-NC dependencies",
    corda2_step2_2_MC_feasibility = "CORDA2 Step 2.2 MC feasibility",
    corda2_step3_HC_OT_dependencies = "CORDA2 Step 3 HC-to-OT dependencies",
    as.character(stage)
  )
}

.rc_corda_stage_context <- function() {
  task <- .rc_layer2_progress_state$current_task
  if (is.list(task) && is.list(task$context)) return(task$context)
  .rc_layer2_task_context(route = "corda2")
}

.rc_corda_stage_progress_dir <- function(stage) {
  task <- .rc_layer2_progress_state$current_task
  base <- if (is.list(task) && !is.null(task$parts_dir)) {
    task$parts_dir
  } else {
    .rc_layer2_progress_dir_from_cache()
  }
  token <- .rc_safe_cache_token(paste(
    .rc_corda_stage_context()$cell_type,
    .rc_corda_stage_context()$medium_scenario,
    stage, Sys.getpid(), sep = "::"
  ))
  file.path(base, paste0("corda_stage_", token))
}

.rc_corda_stage_mark_progress <- function(
    progress_dir, index, total, stage, context, target) {
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)
  marker <- file.path(
    progress_dir, sprintf("done_%09d_%d", as.integer(index), Sys.getpid())
  )
  file.create(marker)
  completed <- min(
    as.integer(total),
    as.integer(length(list.files(progress_dir, pattern = "^done_")))
  )
  remaining <- max(0L, as.integer(total) - completed)
  interval <- max(1L, floor(as.integer(total) / 100L))
  emit <- .rc_progress_enabled(getOption("RegCompassR.progress", TRUE)) &&
    (completed == as.integer(total) || completed <= 2L ||
     (completed %% interval) == 0L)
  .rc_layer2_task_event(
    context = context,
    phase = stage,
    step = completed,
    total = total,
    detail = paste0(
      "completed=", completed,
      "; remaining=", remaining,
      "; current_target=", as.character(target)
    ),
    scope = "corda2_stage",
    run_kind = "primary",
    status = if (remaining == 0L) "complete" else "running",
    parts_dir = dirname(progress_dir),
    emit = emit
  )
  invisible(NULL)
}

.rc_corda_stage_run <- function(targets, stage, FUN) {
  targets <- as.character(targets)
  n_targets <- length(targets)
  label <- .rc_corda_stage_label(stage)
  context <- .rc_corda_stage_context()
  show_progress <- .rc_progress_enabled(getOption("RegCompassR.progress", TRUE))
  if (!n_targets) {
    if (show_progress) {
      message(label, ": no candidate directional targets; skipping")
    }
    .rc_layer2_task_event(
      context, stage, 0L, 1L,
      detail = "no candidate directional targets; remaining=0",
      scope = "corda2_stage", status = "complete", emit = show_progress
    )
    return(list())
  }

  BPPARAM <- .rc_corda_stage_param(n_targets)
  workers <- if (identical(BPPARAM, FALSE)) 1L else
    .rc_layer2_pool_workers(BPPARAM)
  progress_dir <- .rc_corda_stage_progress_dir(stage)
  unlink(progress_dir, recursive = TRUE, force = TRUE)
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit({
    if (!identical(BPPARAM, FALSE) && !is.null(BPPARAM) &&
        requireNamespace("BiocParallel", quietly = TRUE) &&
        methods::is(BPPARAM, "BiocParallelParam") &&
        isTRUE(BiocParallel::bpisup(BPPARAM))) {
      .rc_release_bpparam(BPPARAM)
    }
    unlink(progress_dir, recursive = TRUE, force = TRUE)
    invisible(gc(verbose = FALSE, full = TRUE))
  }, add = TRUE)

  if (show_progress) {
    message(sprintf(
      paste0(
        "%s started | targets=%d | remaining=%d | workers=%d | ",
        "per-worker HiGHS threads=1"
      ),
      label, n_targets, n_targets, workers
    ))
  }
  .rc_layer2_task_event(
    context, stage, 0L, n_targets,
    detail = paste0(
      "stage started; targets=", n_targets,
      "; remaining=", n_targets,
      "; workers=", workers,
      "; solver_threads_per_worker=1"
    ),
    scope = "corda2_stage", status = "running", emit = show_progress
  )

  indexed <- lapply(seq_along(targets), function(i) {
    list(index = as.integer(i), target = targets[[i]])
  })
  worker <- function(item) {
    answer <- FUN(item$target)
    .rc_corda_stage_mark_progress(
      progress_dir = progress_dir,
      index = item$index,
      total = n_targets,
      stage = stage,
      context = context,
      target = item$target
    )
    answer
  }
  answer <- rc_parallel_lapply(indexed, worker, BPPARAM = BPPARAM)

  if (show_progress) {
    message(sprintf(
      "%s complete | completed=%d/%d | remaining=0 | worker pool released",
      label, n_targets, n_targets
    ))
  }
  .rc_layer2_task_event(
    context, stage, n_targets, n_targets,
    detail = paste0(
      "stage complete; completed=", n_targets, "/", n_targets,
      "; remaining=0; worker_pool=released"
    ),
    scope = "corda2_stage", status = "complete", emit = show_progress
  )
  answer
}

.rc_corda2_dependency_target_parallel <- function(
    target, split, directional_class, options, stage, penalized_class,
    solver, time_limit, lower = split$lb, upper = split$ub) {
  engine <- .rc_corda_new_lp_engine(split, solver, time_limit)
  on.exit({ engine <- .rc_corda_release_lp_engine(engine) }, add = TRUE)
  assessed <- .rc_corda2_dependency_assessment_core(
    engine = engine,
    split = split,
    target = target,
    directional_class = directional_class,
    options = options,
    stage = stage,
    penalized_class = penalized_class,
    lower = lower,
    upper = upper
  )
  engine <- assessed$engine
  assessed$engine <- NULL
  assessed$metrics <- .rc_corda_execution_metrics(engine)
  assessed
}

.rc_corda2_maximize_target_parallel <- function(
    target, split, solver, time_limit, lower, upper) {
  engine <- .rc_corda_new_lp_engine(split, solver, time_limit)
  on.exit({ engine <- .rc_corda_release_lp_engine(engine) }, add = TRUE)
  maximum <- .rc_corda2_maximize_target(
    engine, split, target, lower = lower, upper = upper
  )
  engine <- maximum$engine
  maximum$engine <- NULL
  list(
    maximum = maximum,
    metrics = .rc_corda_execution_metrics(engine)
  )
}

.rc_corda2_sum_execution_metrics <- function(parts, n_variables) {
  metrics <- parts[lengths(parts) > 0L]
  sum_field <- function(name) {
    sum(vapply(metrics, function(x) {
      as.numeric(x[[name]] %||% 0)
    }, numeric(1)), na.rm = TRUE)
  }
  all_field <- function(name) {
    length(metrics) > 0L && all(vapply(metrics, function(x) {
      isTRUE(x[[name]])
    }, logical(1)))
  }
  full_values <- sum_field("n_full_vector_numeric_values")
  transmitted <- sum_field("n_transmitted_numeric_values")
  list(
    n_variables = as.integer(n_variables),
    n_solves = as.integer(sum_field("n_solves")),
    n_fallback = as.integer(sum_field("n_fallback")),
    n_objective_coeff_updates = as.integer(
      sum_field("n_objective_coeff_updates")
    ),
    n_bound_index_updates = as.integer(sum_field("n_bound_index_updates")),
    n_sparse_update_calls = as.integer(sum_field("n_sparse_update_calls")),
    n_full_vector_numeric_values = full_values,
    n_transmitted_numeric_values = transmitted,
    n_full_vector_numeric_values_avoided = sum_field(
      "n_full_vector_numeric_values_avoided"
    ),
    transmitted_fraction_of_full = if (full_values > 0) {
      transmitted / full_values
    } else {
      NA_real_
    },
    persistent_solver = all_field("persistent_solver"),
    persistent_disabled = any(vapply(metrics, function(x) {
      isTRUE(x$persistent_disabled)
    }, logical(1))),
    solver_configuration_verified = all_field(
      "solver_configuration_verified"
    ),
    release_policy = "step_local_target_engine_then_worker_pool_release",
    target_parallelism = "within_corda2_stage",
    stage_barrier = TRUE
  )
}

.rc_corda_build_three_stage_core <- function(
    split, classes, options, solver, time_limit) {
  directional_class <- .rc_corda2_directional_class(split, classes)
  initial_directional_class <- directional_class
  hc <- names(directional_class)[directional_class == "HC"]
  mc <- names(directional_class)[directional_class == "MC"]
  nc <- names(directional_class)[directional_class == "NC"]
  ot <- names(directional_class)[directional_class == "OT"]
  inclusion_stage_direction <- stats::setNames(
    rep(NA_character_, length(directional_class)), names(directional_class)
  )
  inclusion_stage_direction[hc] <- "initial_high_confidence"
  task_tables <- list()
  execution_metrics <- list()
  stage_parallelism <- list()

  # Step 1. All targets read the same pre-Step-1 confidence snapshot.
  stage1_hc <- hc
  stage1_mc <- mc
  stage1_nc <- nc
  HCtoMC <- .rc_corda2_empty_dependency_matrix(stage1_hc, stage1_mc)
  HCtoNC <- .rc_corda2_empty_dependency_matrix(stage1_hc, stage1_nc)
  stage1_parts <- .rc_corda_stage_run(
    stage1_hc,
    "corda2_step1_HC_dependencies",
    function(target) {
      .rc_corda2_dependency_target_parallel(
        target = target,
        split = split,
        directional_class = directional_class,
        options = options,
        stage = "corda2_step1_HC_dependencies",
        penalized_class = "stage1",
        solver = solver,
        time_limit = time_limit
      )
    }
  )
  execution_metrics <- c(
    execution_metrics, lapply(stage1_parts, function(x) x$metrics)
  )
  stage_parallelism$step1 <- data.frame(
    stage = "corda2_step1_HC_dependencies",
    n_targets = length(stage1_hc), stringsAsFactors = FALSE
  )
  hc_present <- mc_present <- nc_present <- character()
  blocked_hc <- logical(length(stage1_hc))
  stage1_results <- vector("list", length(stage1_hc))
  for (i in seq_along(stage1_hc)) {
    assessed <- stage1_parts[[i]]
    stage1_results[[i]] <- assessed$result
    if (!isTRUE(assessed$success)) {
      blocked_hc[[i]] <- TRUE
      next
    }
    active <- assessed$active
    hc_present <- union(hc_present, active[directional_class[active] == "HC"])
    mc_used <- assessed$associated[
      directional_class[assessed$associated] == "MC"
    ]
    nc_used <- assessed$associated[
      directional_class[assessed$associated] == "NC"
    ]
    mc_present <- union(mc_present, mc_used)
    nc_present <- union(nc_present, nc_used)
    if (length(mc_used)) HCtoMC[i, mc_used] <- 1L
    if (length(nc_used)) HCtoNC[i, nc_used] <- 1L
  }
  blocked_hc[stage1_hc %in% hc_present] <- FALSE
  retained_stage1_hc <- stage1_hc[!blocked_hc]
  HCtoMC <- HCtoMC[!blocked_hc, , drop = FALSE]
  HCtoNC <- HCtoNC[!blocked_hc, , drop = FALSE]
  promoted_step1_mc <- stage1_mc[stage1_mc %in% mc_present]
  promoted_step1_nc <- stage1_nc[stage1_nc %in% nc_present]
  promoted_step1 <- c(promoted_step1_mc, promoted_step1_nc)
  inclusion_stage_direction[promoted_step1] <-
    "corda2_step1_associated_with_HC"
  directional_class[promoted_step1] <- "HC"
  hc <- c(retained_stage1_hc, promoted_step1)
  mc <- stage1_mc[!stage1_mc %in% mc_present]
  nc <- stage1_nc[!stage1_nc %in% nc_present]
  task_tables$step1 <- .rc_corda2_results_table(stage1_results)
  rm(stage1_parts)
  invisible(gc(verbose = FALSE, full = TRUE))

  # Step 2.1. NC dependency assessment uses one immutable Step-2.1 snapshot.
  stage2_mc_input <- mc
  stage2_nc_input <- nc
  MCxNC <- .rc_corda2_empty_dependency_matrix(stage2_mc_input, stage2_nc_input)
  stage21_parts <- .rc_corda_stage_run(
    stage2_mc_input,
    "corda2_step2_1_MC_NC_dependencies",
    function(target) {
      .rc_corda2_dependency_target_parallel(
        target = target,
        split = split,
        directional_class = directional_class,
        options = options,
        stage = "corda2_step2_1_MC_NC_dependencies",
        penalized_class = "NC",
        solver = solver,
        time_limit = time_limit
      )
    }
  )
  execution_metrics <- c(
    execution_metrics, lapply(stage21_parts, function(x) x$metrics)
  )
  stage_parallelism$step2_1 <- data.frame(
    stage = "corda2_step2_1_MC_NC_dependencies",
    n_targets = length(stage2_mc_input), stringsAsFactors = FALSE
  )
  blocked_mc_step21 <- logical(length(stage2_mc_input))
  stage21_results <- vector("list", length(stage2_mc_input))
  for (i in seq_along(stage2_mc_input)) {
    assessed <- stage21_parts[[i]]
    stage21_results[[i]] <- assessed$result
    if (!isTRUE(assessed$success)) {
      blocked_mc_step21[[i]] <- TRUE
      next
    }
    nc_used <- assessed$associated[
      directional_class[assessed$associated] == "NC"
    ]
    if (length(nc_used)) MCxNC[i, nc_used] <- 1L
  }
  MCxNC <- MCxNC[!blocked_mc_step21, , drop = FALSE]
  mc <- stage2_mc_input[!blocked_mc_step21]
  MCtoNC <- MCxNC
  task_tables$step2_1 <- .rc_corda2_results_table(stage21_results)
  rm(stage21_parts)
  invisible(gc(verbose = FALSE, full = TRUE))

  # Step 2.2. Promotion happens only after all Step-2.1 results are reduced.
  nc_count <- if (ncol(MCxNC)) colSums(MCxNC) else numeric()
  promoted_nc <- names(nc_count)[nc_count >= options$MCxNCthresh]
  if (length(promoted_nc)) {
    promoted_rows <- matrix(
      0L, nrow = length(promoted_nc), ncol = ncol(MCxNC),
      dimnames = list(promoted_nc, colnames(MCxNC))
    )
    MCxNC <- rbind(MCxNC, promoted_rows)
  }
  directional_class[promoted_nc] <- "MC"
  mc <- c(mc, promoted_nc)
  if (length(promoted_nc) && ncol(MCxNC)) {
    MCxNC <- MCxNC[, setdiff(colnames(MCxNC), promoted_nc), drop = FALSE]
  }
  nc <- setdiff(nc, promoted_nc)
  split_step22 <- split
  if (length(nc)) {
    split_step22$lb[nc] <- 0
    split_step22$ub[nc] <- 0
  }
  stage22_parts <- .rc_corda_stage_run(
    mc,
    "corda2_step2_2_MC_feasibility",
    function(target) {
      .rc_corda2_maximize_target_parallel(
        target = target,
        split = split_step22,
        solver = solver,
        time_limit = time_limit,
        lower = split_step22$lb,
        upper = split_step22$ub
      )
    }
  )
  execution_metrics <- c(
    execution_metrics, lapply(stage22_parts, function(x) x$metrics)
  )
  stage_parallelism$step2_2 <- data.frame(
    stage = "corda2_step2_2_MC_feasibility",
    n_targets = length(mc), stringsAsFactors = FALSE
  )
  stage22_results <- vector("list", length(mc))
  blocked_mc_step22 <- logical(length(mc))
  rescue <- vector("list", length(mc))
  for (i in seq_along(mc)) {
    target <- mc[[i]]
    maximum <- stage22_parts[[i]]$maximum
    success <- identical(maximum$answer$status, "optimal") &&
      is.finite(maximum$vmax) && maximum$vmax >= options$flux_threshold
    blocked_mc_step22[[i]] <- !success
    dependencies <- if (target %in% rownames(MCxNC)) {
      colnames(MCxNC)[MCxNC[target, , drop = TRUE] > 0]
    } else {
      character()
    }
    rescue[[i]] <- dependencies
    stage22_results[[i]] <- .rc_corda2_target_result(
      split, target, "corda2_step2_2_MC_feasibility",
      if (success) "optimal" else "target_blocked",
      associated = dependencies,
      target_flux = maximum$vmax,
      vmax = maximum$vmax,
      objective = maximum$answer$objective,
      backend = maximum$answer$backend,
      solver_message = maximum$answer$solver_message %||% "",
      n_solves = 1L,
      opposite = maximum$opposite
    )
  }
  deleted_mc <- mc[blocked_mc_step22]
  rescue_table <- data.frame(
    reaction = deleted_mc,
    dependent_on = vapply(
      rescue[blocked_mc_step22], paste, character(1), collapse = ","
    ),
    stringsAsFactors = FALSE
  )
  feasible_mc <- mc[!blocked_mc_step22]
  inclusion_stage_direction[promoted_nc] <-
    "corda2_step2_2_NC_occurrence_threshold"
  inclusion_stage_direction[feasible_mc[is.na(
    inclusion_stage_direction[feasible_mc]
  )]] <- "corda2_step2_2_MC_feasible"
  directional_class[feasible_mc] <- "HC"
  hc <- c(hc, feasible_mc)
  task_tables$step2_2 <- .rc_corda2_results_table(stage22_results)
  rm(stage22_parts)
  invisible(gc(verbose = FALSE, full = TRUE))

  # Step 3. Remaining MC/NC are blocked before all HC targets are dispatched.
  split_step3 <- split_step22
  allowed_step3 <- union(hc, ot)
  blocked_step3 <- setdiff(colnames(split_step3$S), allowed_step3)
  if (length(blocked_step3)) {
    split_step3$lb[blocked_step3] <- 0
    split_step3$ub[blocked_step3] <- 0
  }
  stage3_parts <- .rc_corda_stage_run(
    hc,
    "corda2_step3_HC_OT_dependencies",
    function(target) {
      .rc_corda2_dependency_target_parallel(
        target = target,
        split = split_step3,
        directional_class = directional_class,
        options = options,
        stage = "corda2_step3_HC_OT_dependencies",
        penalized_class = "OT",
        solver = solver,
        time_limit = time_limit,
        lower = split_step3$lb,
        upper = split_step3$ub
      )
    }
  )
  execution_metrics <- c(
    execution_metrics, lapply(stage3_parts, function(x) x$metrics)
  )
  stage_parallelism$step3 <- data.frame(
    stage = "corda2_step3_HC_OT_dependencies",
    n_targets = length(hc), stringsAsFactors = FALSE
  )
  stage3_results <- vector("list", length(hc))
  ot_present <- character()
  for (i in seq_along(hc)) {
    assessed <- stage3_parts[[i]]
    stage3_results[[i]] <- assessed$result
    if (!isTRUE(assessed$success)) next
    used <- assessed$associated[
      directional_class[assessed$associated] == "OT"
    ]
    ot_present <- union(ot_present, used)
  }
  inclusion_stage_direction[ot_present] <- "corda2_step3_associated_OT"
  included_variables <- unique(c(hc, ot_present))
  task_tables$step3 <- .rc_corda2_results_table(stage3_results)
  rm(stage3_parts)
  invisible(gc(verbose = FALSE, full = TRUE))

  initial_reaction_confidence <- .rc_corda2_reaction_numeric_confidence(
    split, initial_directional_class
  )
  final_reaction_confidence <- .rc_corda2_reaction_numeric_confidence(
    split, initial_directional_class, included_variables
  )
  included_reactions <- unique(as.character(
    split$variable_to_reaction[included_variables]
  ))
  final_reaction_status <- stats::setNames(
    ifelse(names(final_reaction_confidence) %in% included_reactions,
           "included", "excluded"),
    names(final_reaction_confidence)
  )
  inclusion_stage <- stats::setNames(
    rep(NA_character_, length(split$reaction_order)), split$reaction_order
  )
  for (reaction in names(inclusion_stage)) {
    variables <- split$direction_table$variable_id[
      split$direction_table$reaction_id == reaction &
        split$direction_table$variable_id %in% included_variables
    ]
    stages <- inclusion_stage_direction[variables]
    stages <- stages[!is.na(stages)]
    if (length(stages)) inclusion_stage[[reaction]] <- stages[[1L]]
  }

  list(
    included = included_reactions,
    included_directional_variables = included_variables,
    initial_reaction_confidence = initial_reaction_confidence,
    final_reaction_confidence = final_reaction_confidence,
    final_reaction_status = final_reaction_status,
    final_confidence = final_reaction_confidence,
    initial_directional_class = initial_directional_class,
    inclusion_stage = inclusion_stage,
    inclusion_stage_direction = inclusion_stage_direction,
    stage1_associated = unique(as.character(
      split$variable_to_reaction[promoted_step1]
    )),
    stage2_promoted_nc = unique(as.character(
      split$variable_to_reaction[promoted_nc]
    )),
    stage2_promoted_mc = unique(as.character(
      split$variable_to_reaction[feasible_mc]
    )),
    stage3_associated_ot = unique(as.character(
      split$variable_to_reaction[ot_present]
    )),
    blocked_after_stage2 = unique(as.character(
      split$variable_to_reaction[nc]
    )),
    blocked_before_stage3 = unique(as.character(
      split$variable_to_reaction[blocked_step3]
    )),
    blocked_high_confidence_directions = stage1_hc[blocked_hc],
    blocked_medium_directions_step2_1 = stage2_mc_input[blocked_mc_step21],
    blocked_medium_directions_step2_2 = deleted_mc,
    HCtoMC = HCtoMC,
    HCtoNC = HCtoNC,
    MCtoNC = MCtoNC,
    rescue = rescue_table,
    task_diagnostics = .rc_bind_frames_fill(task_tables),
    solver_performance = .rc_corda2_sum_execution_metrics(
      execution_metrics, ncol(split$S)
    ),
    stage_parallelism = .rc_bind_frames_fill(stage_parallelism),
    algorithm = "schultzdre_MATLAB_CORDA2_original_semantics",
    reference_repository = "schultzdre/Constraint-Based-Modeling",
    reference_file = "CORDA2.m",
    stage_update_policy = "original_matlab_directional_order",
    parallel_execution_policy =
      "stage_barrier_parallel_targets_deterministic_ordered_reduce",
    source_semantics = c(
      "split only actively reversible reactions once",
      "close the opposite direction for every tested reaction",
      "fix the target at the original val-or-percentage constraint",
      "step 1 costs MC=sqrt(om), NC=om and all other reactions=1e-3",
      "increase each newly used high-cost reaction once by 1+ci",
      "promote NC directions required by at least MCxNCthresh MC directions",
      "block remaining NC before MC feasibility and add OT only in step 3"
    )
  )
}

# The caller BPPARAM is now a template only. CORDA2 stages own short-lived pools.
.rc_prepare_corda_worker_pool <- function(
    layer2_args, parallel = TRUE, BPPARAM = NULL) {
  state <- list(
    is_corda2 = .rc_layer2_requested_corda2(layer2_args),
    BPPARAM = BPPARAM,
    origin = "not_used",
    started_here = FALSE,
    thread_state = NULL
  )
  if (!isTRUE(state$is_corda2) || !isTRUE(parallel)) return(state)
  if (is.null(state$BPPARAM)) {
    state$BPPARAM <- rc_default_bpparam()
    state$origin <- if (is.null(state$BPPARAM)) {
      "serial_fallback"
    } else {
      "package_default_stage_template"
    }
  } else {
    state$origin <- "caller_supplied_stage_template"
  }
  if (is.null(state$BPPARAM)) return(state)
  if (!requireNamespace("BiocParallel", quietly = TRUE) ||
      !methods::is(state$BPPARAM, "BiocParallelParam")) {
    stop(
      "CORDA2 stage-parallel execution requires a BiocParallelParam object.",
      call. = FALSE
    )
  }
  state
}

.rc_release_corda_worker_pool <- function(state) {
  if (!is.list(state)) return(invisible(NULL))
  if (isTRUE(state$started_here) &&
      requireNamespace("BiocParallel", quietly = TRUE) &&
      methods::is(state$BPPARAM, "BiocParallelParam") &&
      isTRUE(BiocParallel::bpisup(state$BPPARAM))) {
    .rc_release_bpparam(state$BPPARAM)
  }
  if (!is.null(state$thread_state)) {
    .rc_restore_internal_threads(state$thread_state)
  }
  invisible(gc(verbose = FALSE, full = TRUE))
  invisible(NULL)
}

# Move parallelism inward: one model at a time, all workers available per stage.
.rc_corda_should_outer_parallel <- function(n_tasks, pool_workers) FALSE

.rc_build_celltype_medium_corda_cache_serial <-
  .rc_build_celltype_medium_corda_cache
.rc_build_celltype_medium_corda_cache <- function(...) {
  cache <- .rc_build_celltype_medium_corda_cache_serial(...)
  attr(cache, "structural_parallel_task") <-
    "serial_cell_type_x_medium_models_stage_parallel_corda2_targets"
  attr(cache, "structural_parallel_workers") <- 1L
  attr(cache, "structural_dynamic_task_scheduling") <- FALSE
  attr(cache, "corda2_inner_target_parallelism") <- TRUE
  attr(cache, "corda2_stage_barrier_parallelism") <- TRUE
  attr(cache, "corda2_stage_worker_lifecycle") <-
    "start_each_stage_stop_each_stage_full_gc"
  cache
}

.rc_complete_celltype_medium_corda_gem_serial_wrapper <-
  .rc_complete_celltype_medium_corda_gem
.rc_complete_celltype_medium_corda_gem <- function(...) {
  answer <- .rc_complete_celltype_medium_corda_gem_serial_wrapper(...)
  answer$corda_execution$original_matlab$target_parallelism <- TRUE
  answer$corda_execution$original_matlab$parallel_scope <-
    "directional_targets_within_each_original_corda2_stage"
  answer$corda_execution$original_matlab$stage_barrier <- TRUE
  answer$corda_execution$original_matlab$worker_lifecycle <-
    "start_each_stage_stop_each_stage_full_gc"
  answer
}

.rc_layer2_finalize_completion_serial_wrapper <- .rc_layer2_finalize_completion
.rc_layer2_finalize_completion <- function(
    answer, corda_options, is_corda2, solver) {
  answer <- .rc_layer2_finalize_completion_serial_wrapper(
    answer = answer,
    corda_options = corda_options,
    is_corda2 = is_corda2,
    solver = solver
  )
  if (!isTRUE(is_corda2)) return(answer)
  answer$completion_contract$solver_configuration$threads <- 1L
  answer$completion_contract$target_parallelism <-
    "within_each_corda2_stage"
  answer$completion_contract$stage_barrier <- TRUE
  answer$completion_contract$stage_worker_lifecycle <-
    "start_pool_run_all_directional_targets_stop_pool_full_gc_then_next_stage"
  answer$completion_contract$stage_update_policy <-
    "original_matlab_directional_order"
  answer$completion_contract$parallel_execution_policy <-
    "stage_barrier_parallel_targets_deterministic_ordered_reduce"
  answer$params$corda2_inner_target_parallelism <- TRUE
  answer$params$corda2_stage_barrier_parallelism <- TRUE
  answer$union_gem_policy <- paste(
    "one original-CORDA2 reconstruction per cell type and medium;",
    "Step 1, Step 2.1, Step 2.2 and Step 3 remain strict barriers;",
    "directional targets inside each step use the full Layer-2 worker budget"
  )
  answer
}
