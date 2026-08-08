# Stage-barrier parallel helpers for original MATLAB CORDA2 semantics.
#
# These helpers change scheduling only. Every directional target receives the
# same split model, confidence snapshot, objective, bounds and stopping rules as
# the serial implementation. Mathematical state is reduced only at the original
# CORDA2 stage barriers.

.rc_corda_stage_progress_enabled <- function() {
  enabled <- get0(".rc_progress_enabled", mode = "function", inherits = TRUE)
  if (!is.function(enabled)) return(FALSE)
  isTRUE(enabled(getOption("RegCompassR.progress", TRUE)))
}

.rc_corda_stage_event <- function(...) {
  event <- get0(".rc_layer2_task_event", mode = "function", inherits = TRUE)
  if (!is.function(event)) return(invisible(NULL))
  do.call(event, list(...))
}

.rc_corda_stage_algorithm_once <- function(stage) {
  once <- get0(".rc_layer2_algorithm_once", mode = "function", inherits = TRUE)
  if (!is.function(once)) return(invisible(NULL))
  spec <- switch(
    as.character(stage),
    corda2_step1_HC_dependencies = list(
      "corda2_step1", 4L,
      "supporting high-confidence directions with MC/NC dependencies"
    ),
    corda2_step2_1_MC_NC_dependencies = list(
      "corda2_step2_1", 5L,
      "measuring NC dependencies of remaining MC directions"
    ),
    corda2_step2_2_MC_feasibility = list(
      "corda2_step2_2", 6L,
      "testing retained MC directions after NC promotion and blocking"
    ),
    corda2_step3_HC_OT_dependencies = list(
      "corda2_step3", 7L,
      "adding only OT reactions required by retained HC flux"
    ),
    NULL
  )
  if (is.null(spec)) return(invisible(NULL))
  once(spec[[1L]], as.character(stage), spec[[2L]], spec[[3L]])
  invisible(NULL)
}

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
    release <- get0(".rc_release_bpparam", mode = "function", inherits = TRUE)
    if (is.function(release)) {
      release(BPPARAM)
    } else {
      try(BiocParallel::bpstop(BPPARAM), silent = TRUE)
    }
    invisible(gc(verbose = FALSE, full = TRUE))
  }
  invisible(NULL)
}

.rc_corda_stage_param <- function(n_targets) {
  n_targets <- max(0L, as.integer(n_targets[[1L]]))
  context <- get0(
    ".rc_layer2_parallel_context", mode = "environment", inherits = TRUE
  )
  if (n_targets <= 1L || !is.environment(context) ||
      !isTRUE(context$active) || !isTRUE(context$parallel)) return(FALSE)
  template <- context$BPPARAM
  pool_workers <- get0(
    ".rc_corda_pool_workers", mode = "function", inherits = TRUE
  )
  available_workers <- get0(
    "rc_available_workers", mode = "function", inherits = TRUE
  )
  requested <- if (identical(template, FALSE)) {
    1L
  } else if (is.null(template)) {
    if (is.function(available_workers)) {
      max(1L, as.integer(available_workers(default = 1L)))
    } else {
      1L
    }
  } else if (is.function(pool_workers)) {
    pool_workers(template)
  } else if (requireNamespace("BiocParallel", quietly = TRUE) &&
             methods::is(template, "BiocParallelParam")) {
    max(1L, as.integer(BiocParallel::bpnworkers(template)))
  } else {
    1L
  }
  .rc_corda_stage_stop_started_template(template)
  workers <- min(n_targets, max(1L, as.integer(requested)))
  if (workers <= 1L) return(FALSE)
  make_param <- get0("rc_default_bpparam", mode = "function", inherits = TRUE)
  tune_param <- get0(
    ".rc_layer2_tune_task_bpparam", mode = "function", inherits = TRUE
  )
  if (!is.function(make_param)) return(FALSE)
  param <- make_param(
    workers = workers,
    backend = .rc_corda_stage_backend(template)
  )
  if (is.null(param)) return(FALSE)
  if (is.function(tune_param)) param <- tune_param(param, n_targets)
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
  state <- get0(
    ".rc_layer2_progress_state", mode = "environment", inherits = TRUE
  )
  if (is.environment(state)) {
    task <- state$current_task
    if (is.list(task) && is.list(task$context)) return(task$context)
  }
  make_context <- get0(
    ".rc_layer2_task_context", mode = "function", inherits = TRUE
  )
  if (is.function(make_context)) return(make_context(route = "corda2"))
  list(cell_type = "ALL", medium_scenario = "base", route = "corda2")
}

.rc_corda_stage_progress_dir <- function(stage) {
  state <- get0(
    ".rc_layer2_progress_state", mode = "environment", inherits = TRUE
  )
  task <- if (is.environment(state)) state$current_task else NULL
  base <- if (is.list(task) && !is.null(task$parts_dir)) {
    task$parts_dir
  } else {
    progress_dir <- get0(
      ".rc_layer2_progress_dir_from_cache", mode = "function", inherits = TRUE
    )
    if (is.function(progress_dir)) progress_dir() else tempdir()
  }
  context <- .rc_corda_stage_context()
  raw_token <- paste(
    context$cell_type, context$medium_scenario,
    stage, Sys.getpid(), sep = "::"
  )
  safe_token <- get0(".rc_safe_cache_token", mode = "function", inherits = TRUE)
  token <- if (is.function(safe_token)) {
    safe_token(raw_token)
  } else {
    gsub("[^A-Za-z0-9]+", "_", raw_token)
  }
  file.path(base, paste0("corda_stage_", token))
}

.rc_corda_stage_progress_bar <- function(completed, total, width = 28L) {
  total <- max(1L, as.integer(total))
  completed <- max(0L, min(as.integer(completed), total))
  width <- max(10L, as.integer(width))
  filled <- min(width, floor(width * completed / total))
  paste0(
    paste(rep("=", filled), collapse = ""),
    if (filled < width) ">" else "",
    paste(rep(" ", max(0L, width - filled - as.integer(filled < width))),
          collapse = "")
  )
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
  emit <- .rc_corda_stage_progress_enabled() &&
    (completed == as.integer(total) || completed <= 2L ||
     (completed %% interval) == 0L)
  if (emit) {
    message(sprintf(
      "%s [%s] %d/%d (%5.1f%%) remaining=%d current=%s",
      .rc_corda_stage_label(stage),
      .rc_corda_stage_progress_bar(completed, total),
      completed, as.integer(total), 100 * completed / as.integer(total),
      remaining, as.character(target)
    ))
  }
  .rc_corda_stage_event(
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
    emit = FALSE
  )
  invisible(NULL)
}

.rc_corda_stage_chunks <- function(n_targets, workers) {
  n_targets <- max(1L, as.integer(n_targets))
  workers <- max(1L, min(as.integer(workers), n_targets))
  n_chunks <- min(n_targets, max(workers, 16L * workers))
  split(seq_len(n_targets), rep(seq_len(n_chunks), length.out = n_targets))
}

.rc_corda_stage_worker_environment <- function() {
  source_environment <- environment(.rc_corda_stage_worker_environment)
  worker_environment <- new.env(parent = baseenv())
  available <- ls(source_environment, all.names = TRUE)
  keep <- grep("^\\.rc_corda", available, value = TRUE)
  keep <- union(
    keep,
    intersect(
      c(
        "%||%", ".rc_lp_status", ".rc_as_dgCMatrix",
        ".rc_progress_enabled", "rc_validate_gem", "rc_solve_lp"
      ),
      available
    )
  )
  for (name in keep) {
    assign(name, get(name, envir = source_environment, inherits = FALSE),
           envir = worker_environment)
  }
  for (name in keep) {
    value <- get(name, envir = worker_environment, inherits = FALSE)
    if (is.function(value)) {
      environment(value) <- worker_environment
      assign(name, value, envir = worker_environment)
    }
  }
  worker_environment
}

.rc_corda_stage_run <- function(targets, stage, FUN) {
  targets <- as.character(targets)
  n_targets <- length(targets)
  label <- .rc_corda_stage_label(stage)
  context <- .rc_corda_stage_context()
  show_progress <- .rc_corda_stage_progress_enabled()
  if (!n_targets) {
    if (show_progress) {
      message(label, ": no candidate directional targets; skipping")
    }
    .rc_corda_stage_event(
      context, stage, 0L, 1L,
      detail = "no candidate directional targets; remaining=0",
      scope = "corda2_stage", status = "complete", emit = FALSE
    )
    return(list(results = list(), metrics = list(), workers = 1L, chunks = 0L))
  }

  .rc_corda_stage_algorithm_once(stage)
  BPPARAM <- .rc_corda_stage_param(n_targets)
  pool_workers <- get0(
    ".rc_layer2_pool_workers", mode = "function", inherits = TRUE
  )
  workers <- if (identical(BPPARAM, FALSE)) {
    1L
  } else if (is.function(pool_workers)) {
    pool_workers(BPPARAM)
  } else if (requireNamespace("BiocParallel", quietly = TRUE) &&
             methods::is(BPPARAM, "BiocParallelParam")) {
    max(1L, as.integer(BiocParallel::bpnworkers(BPPARAM)))
  } else {
    1L
  }
  chunks <- .rc_corda_stage_chunks(n_targets, workers)
  progress_dir <- .rc_corda_stage_progress_dir(stage)
  unlink(progress_dir, recursive = TRUE, force = TRUE)
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit({
    .rc_corda_stage_stop_started_template(BPPARAM)
    unlink(progress_dir, recursive = TRUE, force = TRUE)
    invisible(gc(verbose = FALSE, full = TRUE))
  }, add = TRUE)

  if (show_progress) {
    message(sprintf(
      paste0(
        "%s started | targets=%d | remaining=%d | workers=%d | chunks=%d | ",
        "per-worker HiGHS threads=1"
      ),
      label, n_targets, n_targets, workers, length(chunks)
    ))
  }
  .rc_corda_stage_event(
    context, stage, 0L, n_targets,
    detail = paste0(
      "stage started; targets=", n_targets,
      "; remaining=", n_targets,
      "; workers=", workers,
      "; chunks=", length(chunks),
      "; solver_threads_per_worker=1"
    ),
    scope = "corda2_stage", status = "running", emit = FALSE
  )

  # Re-home the anonymous stage closure into a plain serializable environment
  # whose parent contains the CORDA2 worker code. This is required for SOCK
  # workers, where nested package closures do not reliably retain indirect
  # internal symbols. Only transport changes; the captured stage snapshot and
  # CORDA2 functions themselves are unchanged.
  worker_namespace <- .rc_corda_stage_worker_environment()
  original_fun_environment <- environment(FUN)
  fun_environment <- new.env(parent = worker_namespace)
  captured_names <- setdiff(
    all.vars(body(FUN), functions = FALSE),
    names(formals(FUN))
  )
  for (name in captured_names) {
    if (exists(name, envir = original_fun_environment, inherits = FALSE)) {
      assign(
        name,
        get(name, envir = original_fun_environment, inherits = FALSE),
        envir = fun_environment
      )
    }
  }
  environment(FUN) <- fun_environment

  worker_environment <- new.env(parent = worker_namespace)
  worker_environment$FUN <- FUN
  worker_environment$targets <- targets
  worker_environment$progress_dir <- progress_dir
  worker_environment$n_targets <- n_targets
  worker_environment$stage <- stage
  worker_environment$context <- context
  worker <- function(index) {
    mark_done <- function(position, target) {
      .rc_corda_stage_mark_progress(
        progress_dir = progress_dir,
        index = position,
        total = n_targets,
        stage = stage,
        context = context,
        target = target
      )
    }
    value <- FUN(targets[index], as.integer(index), mark_done)
    list(index = as.integer(index), value = value)
  }
  environment(worker) <- worker_environment

  if (identical(BPPARAM, FALSE)) {
    parts <- lapply(chunks, worker)
  } else {
    parallel_lapply <- get0(
      "rc_parallel_lapply", mode = "function", inherits = TRUE
    )
    if (!is.function(parallel_lapply)) {
      stop("CORDA2 stage parallelism requires `rc_parallel_lapply`.",
           call. = FALSE)
    }
    parts <- parallel_lapply(chunks, worker, BPPARAM = BPPARAM)
  }

  results <- vector("list", n_targets)
  metrics <- vector("list", length(parts))
  for (i in seq_along(parts)) {
    part <- parts[[i]]
    value <- part$value
    if (!is.list(value) || !is.list(value$results) ||
        length(value$results) != length(part$index)) {
      stop("A CORDA2 stage worker returned malformed target results.",
           call. = FALSE)
    }
    results[part$index] <- value$results
    metrics[[i]] <- value$metrics %||% list()
  }

  .rc_corda_stage_stop_started_template(BPPARAM)
  if (show_progress) {
    message(sprintf(
      "%s complete | completed=%d/%d | remaining=0 | worker pool released",
      label, n_targets, n_targets
    ))
  }
  .rc_corda_stage_event(
    context, stage, n_targets, n_targets,
    detail = paste0(
      "stage complete; completed=", n_targets, "/", n_targets,
      "; remaining=0; worker_pool=released"
    ),
    scope = "corda2_stage", status = "complete", emit = FALSE
  )
  list(
    results = results,
    metrics = metrics,
    workers = as.integer(workers),
    chunks = as.integer(length(chunks))
  )
}

.rc_corda2_dependency_chunk_parallel <- function(
    targets, indices, mark_done, split, directional_class, options,
    stage, penalized_class, solver, time_limit,
    lower = split$lb, upper = split$ub) {
  engine <- .rc_corda_new_lp_engine(split, solver, time_limit)
  on.exit({ engine <- .rc_corda_release_lp_engine(engine) }, add = TRUE)
  results <- vector("list", length(targets))
  for (i in seq_along(targets)) {
    assessed <- .rc_corda2_dependency_assessment_core(
      engine = engine,
      split = split,
      target = targets[[i]],
      directional_class = directional_class,
      options = options,
      stage = stage,
      penalized_class = penalized_class,
      lower = lower,
      upper = upper
    )
    engine <- assessed$engine
    assessed$engine <- NULL
    results[[i]] <- assessed
    mark_done(indices[[i]], targets[[i]])
  }
  list(
    results = results,
    metrics = .rc_corda_execution_metrics(engine)
  )
}

.rc_corda2_maximize_chunk_parallel <- function(
    targets, indices, mark_done, split, solver, time_limit, lower, upper) {
  engine <- .rc_corda_new_lp_engine(split, solver, time_limit)
  on.exit({ engine <- .rc_corda_release_lp_engine(engine) }, add = TRUE)
  results <- vector("list", length(targets))
  for (i in seq_along(targets)) {
    maximum <- .rc_corda2_maximize_target(
      engine, split, targets[[i]], lower = lower, upper = upper
    )
    engine <- maximum$engine
    maximum$engine <- NULL
    results[[i]] <- maximum
    mark_done(indices[[i]], targets[[i]])
  }
  list(
    results = results,
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
    release_policy = "stage_chunk_engine_then_worker_pool_release",
    target_parallelism = "within_corda2_stage",
    stage_barrier = TRUE
  )
}
