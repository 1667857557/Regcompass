# Keep one worker pool alive across independent CORDA2 model instances.

.rc_layer2_requested_corda2 <- function(layer2_args) {
  model_params <- if (is.list(layer2_args)) {
    layer2_args$model_params %||% list()
  } else {
    list()
  }
  requested <- as.character(model_params$model_completion %||% "corda2")
  length(requested) == 1L && !is.na(requested) &&
    requested %in% c("corda2", "corda")
}

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
      "package_default"
    }
  } else {
    state$origin <- "caller_supplied"
  }
  if (is.null(state$BPPARAM)) return(state)
  if (!requireNamespace("BiocParallel", quietly = TRUE) ||
      !methods::is(state$BPPARAM, "BiocParallelParam")) {
    stop(
      "CORDA2 outer-model parallel execution requires a ",
      "BiocParallelParam object.",
      call. = FALSE
    )
  }
  if (!isTRUE(BiocParallel::bpisup(state$BPPARAM))) {
    state$thread_state <- .rc_set_internal_single_thread()
    start_error <- tryCatch({
      BiocParallel::bpstart(state$BPPARAM)
      NULL
    }, error = function(e) e)
    if (inherits(start_error, "error")) {
      .rc_restore_internal_threads(state$thread_state)
      stop(
        "Unable to start the CORDA2 BiocParallel worker pool: ",
        conditionMessage(start_error),
        call. = FALSE
      )
    }
    state$started_here <- TRUE
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

# Hierarchical Layer 2 progress reporting -------------------------------------
#
# Layer 2 has two distinct progress scopes:
#   1. one controller-owned 12-part stage bar for the complete workflow;
#   2. task-scoped events for each cell-type-by-medium structural model and
#      scoring pass. Worker events are written to per-process files and merged
#      into layer2_task_progress.tsv when the stage finishes or fails.
#
# Canonical Layer 2 functions call these shared progress helpers directly.
# No runtime function capture or same-name redefinition is used.

.rc_layer2_progress_parts <- 12L
.rc_layer2_progress_state <- new.env(parent = emptyenv())
.rc_layer2_progress_state$active <- FALSE
.rc_layer2_progress_state$owner_pid <- NA_integer_
.rc_layer2_progress_state$monitor <- NULL
.rc_layer2_progress_state$outdir <- NULL
.rc_layer2_progress_state$parts_dir <- NULL
.rc_layer2_progress_state$run_kind <- "primary"
.rc_layer2_progress_state$task_sequence <- 0L
.rc_layer2_progress_state$task_started <- new.env(parent = emptyenv())
.rc_layer2_progress_state$task_step <- new.env(parent = emptyenv())
.rc_layer2_progress_state$current_task <- NULL
.rc_layer2_progress_state$inside_dependency <- FALSE
.rc_layer2_progress_state$algorithm_flags <- new.env(parent = emptyenv())


.rc_layer2_progress_sanitize <- function(x, empty = "NA") {
  value <- if (is.null(x) || !length(x)) empty else {
    paste(as.character(x), collapse = ",")
  }
  value <- gsub("[\t\r\n;]+", " ", value)
  value <- trimws(value)
  if (!nzchar(value)) empty else value
}

.rc_layer2_progress_dir_from_cache <- function(cache_dir = NULL) {
  if (!is.null(cache_dir) && length(cache_dir)) {
    value <- normalizePath(
      as.character(cache_dir[[1L]]), winslash = "/", mustWork = FALSE
    )
    if (identical(basename(value), "corda2")) value <- dirname(value)
    parent <- dirname(dirname(value))
    if (nzchar(parent) && !identical(parent, ".")) {
      return(file.path(parent, "layer2_task_progress_parts"))
    }
  }
  state <- .rc_layer2_progress_state
  if (isTRUE(state$active) && !is.null(state$parts_dir)) {
    return(state$parts_dir)
  }
  file.path(tempdir(), "RegCompassR_layer2_task_progress_parts")
}

.rc_layer2_progress_cache_dir_from_frames <- function() {
  frames <- rev(sys.frames())
  for (frame in frames) {
    if (!exists("cache_dir", envir = frame, inherits = FALSE)) next
    value <- tryCatch(
      get("cache_dir", envir = frame, inherits = FALSE),
      error = function(error) NULL
    )
    if (!is.null(value) && length(value)) return(as.character(value[[1L]]))
  }
  NULL
}

.rc_layer2_progress_begin <- function(monitor) {
  if (is.null(monitor) || !is.environment(monitor)) return(invisible(NULL))
  state <- .rc_layer2_progress_state
  state$active <- TRUE
  state$owner_pid <- Sys.getpid()
  state$monitor <- monitor
  state$outdir <- monitor$outdir
  state$parts_dir <- if (is.null(monitor$outdir)) {
    file.path(tempdir(), "RegCompassR_layer2_task_progress_parts")
  } else {
    file.path(monitor$outdir, "layer2_task_progress_parts")
  }
  dir.create(state$parts_dir, recursive = TRUE, showWarnings = FALSE)
  old_parts <- list.files(state$parts_dir, full.names = TRUE, all.files = TRUE)
  old_parts <- old_parts[!basename(old_parts) %in% c(".", "..")]
  if (length(old_parts)) unlink(old_parts, recursive = TRUE, force = TRUE)
  state$run_kind <- "primary"
  state$task_sequence <- 0L
  state$task_started <- new.env(parent = emptyenv())
  state$task_step <- new.env(parent = emptyenv())
  state$current_task <- NULL
  state$inside_dependency <- FALSE
  state$algorithm_flags <- new.env(parent = emptyenv())
  invisible(state)
}

.rc_layer2_progress_reset <- function() {
  state <- .rc_layer2_progress_state
  state$active <- FALSE
  state$owner_pid <- NA_integer_
  state$monitor <- NULL
  state$outdir <- NULL
  state$parts_dir <- NULL
  state$run_kind <- "primary"
  state$current_task <- NULL
  state$inside_dependency <- FALSE
  state$task_started <- new.env(parent = emptyenv())
  state$task_step <- new.env(parent = emptyenv())
  invisible(NULL)
}

.rc_layer2_progress_is_controller <- function() {
  state <- .rc_layer2_progress_state
  isTRUE(state$active) && identical(Sys.getpid(), state$owner_pid) &&
    is.environment(state$monitor)
}

.rc_layer2_overall_event <- function(
    phase, current, detail = NULL, context = list(), status = "running") {
  if (!.rc_layer2_progress_is_controller()) return(invisible(NULL))
  .rc_step_monitor_event(
    .rc_layer2_progress_state$monitor,
    phase = phase,
    detail = detail,
    current = current,
    total = .rc_layer2_progress_parts,
    context = context,
    status = status,
    emit = TRUE
  )
}

.rc_layer2_medium_id <- function(medium_table) {
  if (is.null(medium_table) || !is.data.frame(medium_table) ||
      !nrow(medium_table)) return("base")
  column <- if ("medium_scenario_id" %in% colnames(medium_table)) {
    medium_table$medium_scenario_id
  } else {
    "custom"
  }
  values <- unique(trimws(as.character(column)))
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values)) paste(values, collapse = ",") else "base"
}

.rc_layer2_task_key <- function(context, scope, run_kind) {
  paste(
    .rc_layer2_progress_sanitize(context$cell_type, "ALL"),
    .rc_layer2_progress_sanitize(context$medium_scenario, "base"),
    .rc_layer2_progress_sanitize(context$route, "unknown"),
    .rc_layer2_progress_sanitize(scope, "task"),
    .rc_layer2_progress_sanitize(run_kind, "primary"),
    sep = "\001"
  )
}

.rc_layer2_task_context <- function(
    cell_type = "ALL", medium_scenario = "base", route = "unknown") {
  list(
    cell_type = .rc_layer2_progress_sanitize(cell_type, "ALL"),
    medium_scenario = .rc_layer2_progress_sanitize(
      medium_scenario, "base"
    ),
    route = .rc_layer2_progress_sanitize(route, "unknown")
  )
}

.rc_layer2_task_progress_file <- function(parts_dir) {
  dir.create(parts_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(parts_dir, paste0("events_", Sys.getpid(), ".tsv"))
}

.rc_layer2_task_event <- function(
    context, phase, step, total, detail = NULL,
    scope = "structural", run_kind = NULL,
    status = "running", parts_dir = NULL, emit = NULL) {
  if (!is.list(context)) {
    stop("Layer 2 task progress context must be a list.", call. = FALSE)
  }
  step <- suppressWarnings(as.integer(step[[1L]]))
  total <- suppressWarnings(as.integer(total[[1L]]))
  if (!is.finite(total) || total < 1L) total <- 1L
  if (!is.finite(step)) step <- 0L
  step <- max(0L, min(step, total))
  state <- .rc_layer2_progress_state
  run_kind <- .rc_layer2_progress_sanitize(
    run_kind %||% state$run_kind %||% "primary", "primary"
  )
  scope <- .rc_layer2_progress_sanitize(scope, "task")
  parts_dir <- parts_dir %||% state$parts_dir %||%
    .rc_layer2_progress_dir_from_cache()
  dir.create(parts_dir, recursive = TRUE, showWarnings = FALSE)
  key <- .rc_layer2_task_key(context, scope, run_kind)
  if (!exists(key, envir = state$task_started, inherits = FALSE)) {
    assign(key, unname(proc.time()[["elapsed"]]),
           envir = state$task_started)
  }
  started <- get(key, envir = state$task_started, inherits = FALSE)
  elapsed <- max(0, unname(proc.time()[["elapsed"]]) - as.numeric(started))
  state$task_sequence <- as.integer(state$task_sequence %||% 0L) + 1L
  assign(key, step, envir = state$task_step)
  row <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3%z"),
    pid = Sys.getpid(),
    sequence = state$task_sequence,
    run_kind = run_kind,
    scope = scope,
    cell_type = .rc_layer2_progress_sanitize(context$cell_type, "ALL"),
    medium_scenario = .rc_layer2_progress_sanitize(
      context$medium_scenario, "base"
    ),
    route = .rc_layer2_progress_sanitize(context$route, "unknown"),
    phase = .rc_layer2_progress_sanitize(phase, "unknown"),
    status = .rc_layer2_progress_sanitize(status, "running"),
    step = step,
    total = total,
    percent = 100 * step / total,
    elapsed_seconds = elapsed,
    elapsed_hms = .rc_format_elapsed(elapsed),
    detail = .rc_layer2_progress_sanitize(detail, ""),
    stringsAsFactors = FALSE
  )
  path <- .rc_layer2_task_progress_file(parts_dir)
  append <- file.exists(path)
  utils::write.table(
    row, path, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = !append, append = append
  )
  if (is.null(emit)) {
    emit <- .rc_progress_enabled(getOption("RegCompassR.progress", TRUE))
  }
  if (isTRUE(emit)) {
    message(sprintf(
      paste0(
        "RegCompass layer2 task ",
        "[cell_type=%s | medium=%s | route=%s | run=%s | scope=%s] ",
        "%d/%d (%5.1f%%) elapsed=%s phase=%s%s"
      ),
      row$cell_type, row$medium_scenario, row$route, row$run_kind,
      row$scope, row$step, row$total, row$percent, row$elapsed_hms,
      row$phase,
      if (nzchar(row$detail)) paste0(" | ", row$detail) else ""
    ))
  }
  invisible(row)
}

.rc_layer2_task_last_step <- function(
    context, scope = "structural", run_kind = "primary", default = 1L) {
  key <- .rc_layer2_task_key(context, scope, run_kind)
  steps <- .rc_layer2_progress_state$task_step
  if (exists(key, envir = steps, inherits = FALSE)) {
    return(as.integer(get(key, envir = steps, inherits = FALSE)))
  }
  as.integer(default)
}

.rc_layer2_task_push <- function(context, route, total, parts_dir = NULL) {
  state <- .rc_layer2_progress_state
  previous <- list(
    current_task = state$current_task,
    inside_dependency = state$inside_dependency,
    algorithm_flags = state$algorithm_flags
  )
  resolved_parts_dir <- parts_dir
  if (is.null(resolved_parts_dir)) {
    resolved_cache_dir <- .rc_layer2_progress_cache_dir_from_frames()
    resolved_parts_dir <- .rc_layer2_progress_dir_from_cache(
      resolved_cache_dir
    )
  }
  state$current_task <- list(
    context = context,
    route = route,
    total = as.integer(total),
    parts_dir = resolved_parts_dir
  )
  state$inside_dependency <- FALSE
  state$algorithm_flags <- new.env(parent = emptyenv())
  .rc_layer2_task_event(
    context, "task_start", 1L, total,
    detail = "task inputs selected",
    scope = "structural", run_kind = "primary",
    parts_dir = state$current_task$parts_dir
  )
  previous
}

.rc_layer2_task_pop <- function(previous) {
  state <- .rc_layer2_progress_state
  state$current_task <- previous$current_task
  state$inside_dependency <- previous$inside_dependency
  state$algorithm_flags <- previous$algorithm_flags
  invisible(NULL)
}

.rc_layer2_current_task_event <- function(
    phase, step, detail = NULL, status = "running") {
  task <- .rc_layer2_progress_state$current_task
  if (is.null(task) || !is.list(task)) return(invisible(NULL))
  .rc_layer2_task_event(
    task$context, phase, step, task$total, detail = detail,
    scope = "structural", run_kind = "primary", status = status,
    parts_dir = task$parts_dir
  )
}

.rc_layer2_algorithm_once <- function(flag, phase, step, detail = NULL) {
  flags <- .rc_layer2_progress_state$algorithm_flags
  if (exists(flag, envir = flags, inherits = FALSE)) return(invisible(NULL))
  assign(flag, TRUE, envir = flags)
  .rc_layer2_current_task_event(phase, step, detail)
}

.rc_layer2_cache_progress_dir <- function(model_cache) {
  if (is.list(model_cache) && length(model_cache)) {
    entry <- model_cache[[1L]]
    file <- as.character(entry$file %||% "")
    if (nzchar(file)) {
      directory <- dirname(file)
      if (identical(basename(directory), "corda2")) {
        directory <- dirname(directory)
      }
      return(.rc_layer2_progress_dir_from_cache(directory))
    }
  }
  .rc_layer2_progress_dir_from_cache()
}

.rc_layer2_model_contexts <- function(model_cache, mode = NULL) {
  if (!is.list(model_cache) || !length(model_cache)) return(list())
  files <- vapply(model_cache, function(entry) {
    as.character(entry$file %||% "memory")
  }, character(1))
  representatives <- match(unique(files), files)
  lapply(representatives, function(index) {
    entry <- model_cache[[index]]
    strategy <- as.character(entry$build_strategy %||% "")
    route <- if (identical(mode, "full_gem")) {
      "full_gem"
    } else if (grepl("corda2", strategy, ignore.case = TRUE)) {
      "corda2"
    } else {
      "fastcore"
    }
    context <- .rc_layer2_task_context(
      cell_type = entry$cell_type %||% "ALL",
      medium_scenario = entry$medium_scenario %||% "base",
      route = route
    )
    list(
      context = context,
      file = files[[index]],
      n_targets = sum(files == files[[index]])
    )
  })
}

.rc_layer2_collect_task_progress <- function(outdir) {
  if (is.null(outdir) || !length(outdir)) return(invisible(NULL))
  parts_dir <- file.path(outdir, "layer2_task_progress_parts")
  files <- list.files(
    parts_dir, pattern = "^events_[0-9]+\\.tsv$", full.names = TRUE
  )
  if (!length(files)) return(invisible(NULL))
  rows <- lapply(files, function(file) {
    tryCatch(
      utils::read.delim(file, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(error) data.frame()
    )
  })
  rows <- rows[vapply(rows, nrow, integer(1)) > 0L]
  if (!length(rows)) return(invisible(NULL))
  events <- .rc_bind_frames_fill(rows)
  order_columns <- intersect(
    c("timestamp", "pid", "sequence"), colnames(events)
  )
  if (length(order_columns)) {
    events <- events[do.call(order, events[order_columns]), , drop = FALSE]
  }
  rownames(events) <- NULL
  utils::write.table(
    events,
    file.path(outdir, "layer2_task_progress.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
  )
  key_columns <- intersect(
    c("run_kind", "scope", "cell_type", "medium_scenario", "route"),
    colnames(events)
  )
  if (length(key_columns)) {
    key <- do.call(paste, c(events[key_columns], sep = "\001"))
    latest <- !duplicated(key, fromLast = TRUE)
    status <- events[latest, , drop = FALSE]
    status <- status[do.call(order, status[key_columns]), , drop = FALSE]
    utils::write.table(
      status,
      file.path(outdir, "layer2_task_status.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
    )
  }
  invisible(events)
}

# Controller monitor: Layer 2 starts with a meaningful total rather than 0/1.
# Overall primary and RNA-control phases.
# Structural FASTCORE task phases.
# Original CORDA2 task phases.
# Full-GEM model construction has no cell-type partition.
# Per-model vmax and Step 2 scoring phases.
# Structural cache completion advances the controller-owned overall bar.
# Final comparison and validation phases.
