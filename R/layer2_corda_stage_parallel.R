# Stage-barrier parallel helpers for original MATLAB CORDA2 semantics.
#
# These helpers change scheduling only. Every directional target receives the
# same split model, confidence snapshot, objective, bounds and stopping rules as
# the serial implementation. Mathematical state is reduced only at the original
# CORDA2 stage barriers.

.rc_corda_stage_backend <- function(BPPARAM) {
  if (identical(BPPARAM, FALSE) || is.null(BPPARAM) ||
      !requireNamespace("BiocParallel", quietly = TRUE) ||
      !methods::is(BPPARAM, "BiocParallelParam")) return("auto")
  if (methods::is(BPPARAM, "SnowParam")) return("snow")
  if (methods::is(BPPARAM, "MulticoreParam")) return("multicore")
  "auto"
}

.rc_corda_stage_stop_started_template <- function(BPPARAM) {
  if (identical(BPPARAM, FALSE) || is.null(BPPARAM) ||
      !requireNamespace("BiocParallel", quietly = TRUE) ||
      !methods::is(BPPARAM, "BiocParallelParam")) return(invisible(NULL))
  if (isTRUE(BiocParallel::bpisup(BPPARAM))) {
    .rc_release_bpparam(BPPARAM)
    invisible(gc(verbose = FALSE, full = TRUE))
  }
  invisible(NULL)
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
  .rc_corda_stage_stop_started_template(template)
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
  context <- .rc_corda_stage_context()
  token <- .rc_safe_cache_token(paste(
    context$cell_type, context$medium_scenario,
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

  .rc_corda_stage_stop_started_template(BPPARAM)
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
