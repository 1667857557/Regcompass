# Keep one worker pool alive across independent CORDA2 model instances.

.rc_layer2_requested_corda2 <- function(layer2_args) {
  model_params <- if (is.list(layer2_args)) {
    layer2_args$model_params %||% list()
  } else {
    list()
  }
  requested <- as.character(model_params$model_completion %||% "fastcore")
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
# This section is loaded after every Layer 2 implementation so the wrappers
# below observe the final runtime functions, including the reaction-parallel
# scoring engine and the original CORDA2 implementation.

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

.rc_layer2_progress_original <- new.env(parent = emptyenv())

.rc_layer2_progress_capture <- function(name) {
  value <- get0(name, mode = "function", inherits = TRUE)
  if (!is.function(value)) return(NULL)
  assign(name, value, envir = .rc_layer2_progress_original)
  value
}

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
    value <- get0("cache_dir", envir = frame, inherits = TRUE)
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
  state$current_task <- list(
    context = context,
    route = route,
    total = as.integer(total),
    parts_dir = parts_dir %||% .rc_layer2_progress_dir_from_cache(
      .rc_layer2_progress_cache_dir_from_frames()
    )
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
.rc_layer2_progress_capture(".rc_step_monitor_start")
.rc_step_monitor_start <- function(
    stage, outdir, progress = TRUE, total_parts = 1L) {
  if (identical(as.character(stage[[1L]]), "layer2") &&
      identical(as.integer(total_parts[[1L]]), 1L)) {
    total_parts <- .rc_layer2_progress_parts
  }
  monitor <- .rc_layer2_progress_original$.rc_step_monitor_start(
    stage, outdir, progress, total_parts
  )
  if (identical(as.character(stage[[1L]]), "layer2")) {
    .rc_layer2_progress_begin(monitor)
  }
  monitor
}

.rc_layer2_progress_capture(".rc_step_monitor_finish")
.rc_step_monitor_finish <- function(
    value, monitor, status = "success", details = NULL) {
  if (!is.null(monitor) && is.environment(monitor) &&
      identical(as.character(monitor$timer$stage), "layer2")) {
    .rc_layer2_collect_task_progress(monitor$outdir)
  }
  .rc_layer2_progress_original$.rc_step_monitor_finish(
    value, monitor, status, details
  )
}

.rc_layer2_progress_capture(".rc_step_monitor_fail")
.rc_step_monitor_fail <- function(monitor) {
  is_layer2 <- !is.null(monitor) && is.environment(monitor) &&
    identical(as.character(monitor$timer$stage), "layer2")
  if (is_layer2) .rc_layer2_collect_task_progress(monitor$outdir)
  answer <- .rc_layer2_progress_original$.rc_step_monitor_fail(monitor)
  if (is_layer2) .rc_layer2_progress_reset()
  answer
}

# Overall primary and RNA-control phases.
.rc_layer2_progress_capture(".rc_run_microcompass_monitored")
.rc_run_microcompass_monitored <- function(..., progress_monitor = NULL) {
  args <- list(...)
  state <- .rc_layer2_progress_state
  previous_run <- state$run_kind
  state$run_kind <- "primary"
  on.exit({ state$run_kind <- previous_run }, add = TRUE)
  media <- unique(as.character(
    args$medium_scenarios$medium_scenario_id %||% "base"
  ))
  celltypes <- if (identical(as.character(args$mode), "meta_module_gem") &&
                   is.data.frame(args$reaction_membership)) {
    column <- as.character(args$celltype_col %||% "cell_type")
    unique(as.character(args$reaction_membership[[column]]))
  } else {
    "ALL"
  }
  completion <- .rc_layer2_completion_context$model_completion %||%
    if (identical(as.character(args$mode), "full_gem")) "none" else "fastcore"
  target_count <- if (is.data.frame(args$target_reactions)) {
    length(unique(as.character(args$target_reactions$reaction_id)))
  } else {
    length(unique(as.character(args$target_reactions)))
  }
  .rc_layer2_overall_event(
    "layer2_plan_ready", 1L,
    detail = paste(
      "execution plan resolved:", length(celltypes), "cell types x",
      length(media), "media;", target_count, "target reactions"
    ),
    context = list(
      model_mode = args$mode,
      completion = completion,
      cell_types = length(celltypes),
      media = length(media),
      parallel = isTRUE(args$parallel %||% TRUE)
    )
  )
  .rc_layer2_overall_event(
    "primary_engine_start", 2L,
    "starting structural construction and primary multiome scoring"
  )
  answer <- do.call(rc_run_microcompass, args)
  .rc_layer2_overall_event(
    "primary_scoring_complete", 5L,
    "primary multiome Step 2 scoring completed"
  )
  .rc_layer2_overall_event(
    "primary_engine_complete", 6L,
    "primary structural models and directional scores assembled"
  )
  answer
}

.rc_layer2_progress_capture(".rc_run_microcompass_engine_monitored")
.rc_run_microcompass_engine_monitored <- function(
    ..., progress_monitor = NULL) {
  args <- list(...)
  state <- .rc_layer2_progress_state
  previous_run <- state$run_kind
  state$run_kind <- "rna_control"
  on.exit({ state$run_kind <- previous_run }, add = TRUE)
  .rc_layer2_overall_event(
    "rna_control_start", 7L,
    "reusing structural models for RNA-only control scoring"
  )
  answer <- do.call(.rc_run_microcompass_engine, args)
  .rc_layer2_overall_event(
    "rna_control_complete", 8L,
    "RNA-only control scoring completed on the shared structural cache"
  )
  answer
}

# Structural FASTCORE task phases.
.rc_layer2_progress_capture(".rc_complete_celltype_medium_union_gem")
.rc_complete_celltype_medium_union_gem <- function(...) {
  args <- list(...)
  context <- .rc_layer2_task_context(
    cell_type = args$cell_type %||% args[[4L]],
    medium_scenario = .rc_layer2_medium_id(
      args$medium_table %||% if (length(args) >= 5L) args[[5L]] else NULL
    ),
    route = "fastcore"
  )
  previous <- .rc_layer2_task_push(context, "fastcore", 6L)
  on.exit(.rc_layer2_task_pop(previous), add = TRUE)
  tryCatch({
    answer <- do.call(
      .rc_layer2_progress_original$.rc_complete_celltype_medium_union_gem,
      args
    )
    if (!exists(
          "fastcore_support",
          envir = .rc_layer2_progress_state$algorithm_flags,
          inherits = FALSE
        )) {
      .rc_layer2_algorithm_once(
        "fastcore_support", "fastcore_support_completion", 5L,
        "skipped because all parent-feasible targets were already supported"
      )
    }
    .rc_layer2_current_task_event(
      "fastcore_model_ready", 6L,
      detail = paste0(
        "reactions=", ncol(answer$S),
        "; support=",
        answer$build_params$n_celltype_fastcore_support_reactions %||% NA_integer_
      ),
      status = "complete"
    )
    answer
  }, error = function(error) {
    .rc_layer2_current_task_event(
      "fastcore_task_error",
      .rc_layer2_task_last_step(context),
      conditionMessage(error), status = "error"
    )
    stop(error)
  })
}

.rc_layer2_progress_capture(".rc_fastcore_parent")
.rc_fastcore_parent <- function(...) {
  task <- .rc_layer2_progress_state$current_task
  if (!is.null(task) && identical(task$route, "fastcore")) {
    .rc_layer2_current_task_event(
      "medium_parent_start", 2L,
      "applying medium bounds and checking parent feasibility"
    )
  }
  answer <- do.call(.rc_layer2_progress_original$.rc_fastcore_parent, list(...))
  if (!is.null(task) && identical(task$route, "fastcore")) {
    .rc_layer2_current_task_event(
      "fastcc_complete", 3L,
      detail = paste0(
        "consistent=", length(answer$fastcc_consistent_reactions %||% character()),
        "; inconsistent=", length(answer$fastcc_inconsistent_reactions %||% character())
      )
    )
  }
  answer
}

.rc_layer2_progress_capture(".rc_directional_feasibility")
.rc_directional_feasibility <- function(...) {
  task <- .rc_layer2_progress_state$current_task
  if (!is.null(task) && identical(task$route, "fastcore")) {
    .rc_layer2_algorithm_once(
      "fastcore_direction_scan", "core_direction_scan", 4L,
      "testing parent and biological-submodel target directions"
    )
  }
  do.call(.rc_layer2_progress_original$.rc_directional_feasibility, list(...))
}

.rc_layer2_progress_capture(".rc_fastcore_complete_direction")
.rc_fastcore_complete_direction <- function(...) {
  args <- list(...)
  task <- .rc_layer2_progress_state$current_task
  if (!is.null(task) && identical(task$route, "fastcore")) {
    direction <- args$direction %||%
      if (length(args) >= 5L) args[[5L]] else "unknown"
    .rc_layer2_algorithm_once(
      "fastcore_support", "fastcore_support_completion", 5L,
      detail = paste0("target_direction=", direction)
    )
  }
  do.call(
    .rc_layer2_progress_original$.rc_fastcore_complete_direction,
    args
  )
}

# Original CORDA2 task phases.
.rc_layer2_progress_capture(".rc_complete_celltype_medium_corda_gem")
.rc_complete_celltype_medium_corda_gem <- function(...) {
  args <- list(...)
  context <- .rc_layer2_task_context(
    cell_type = args$cell_type %||% args[[4L]],
    medium_scenario = .rc_layer2_medium_id(args$medium_table),
    route = "corda2"
  )
  previous <- .rc_layer2_task_push(context, "corda2", 9L)
  on.exit(.rc_layer2_task_pop(previous), add = TRUE)
  tryCatch({
    answer <- do.call(
      .rc_layer2_progress_original$.rc_complete_celltype_medium_corda_gem,
      args
    )
    .rc_layer2_current_task_event(
      "corda2_model_ready", 9L,
      detail = paste0(
        "included_reactions=", ncol(answer$S),
        "; lp_solves=",
        answer$corda_reconstruction$solver_performance$n_solves %||% NA_integer_
      ),
      status = "complete"
    )
    answer
  }, error = function(error) {
    .rc_layer2_current_task_event(
      "corda2_task_error",
      .rc_layer2_task_last_step(context),
      conditionMessage(error), status = "error"
    )
    stop(error)
  })
}

.rc_layer2_progress_capture(".rc_corda_parent")
.rc_corda_parent <- function(...) {
  task <- .rc_layer2_progress_state$current_task
  if (!is.null(task) && identical(task$route, "corda2")) {
    .rc_layer2_current_task_event(
      "medium_parent_start", 2L,
      "applying COMPASS medium bounds to the complete parent GEM"
    )
  }
  answer <- do.call(.rc_layer2_progress_original$.rc_corda_parent, list(...))
  if (!is.null(task) && identical(task$route, "corda2")) {
    .rc_layer2_current_task_event(
      "medium_parent_complete", 2L,
      detail = paste0(
        "reactions=", answer$corda_parent_n_reactions %||% ncol(answer$S),
        "; open=", answer$corda_parent_n_open_reactions %||% NA_integer_
      )
    )
  }
  answer
}

.rc_layer2_progress_capture(".rc_corda_classify_reactions")
.rc_corda_classify_reactions <- function(...) {
  task <- .rc_layer2_progress_state$current_task
  if (!is.null(task) && identical(task$route, "corda2")) {
    .rc_layer2_current_task_event(
      "corda2_confidence_mapping", 3L,
      "mapping multiome evidence to HC, MC, NC and OT"
    )
  }
  do.call(.rc_layer2_progress_original$.rc_corda_classify_reactions, list(...))
}

.rc_layer2_progress_capture(".rc_corda2_dependency_assessment")
.rc_corda2_dependency_assessment <- function(...) {
  args <- list(...)
  stage <- args$stage %||% if (length(args) >= 6L) args[[6L]] else ""
  task <- .rc_layer2_progress_state$current_task
  if (!is.null(task) && identical(task$route, "corda2")) {
    if (identical(stage, "corda2_step1_HC_dependencies")) {
      .rc_layer2_algorithm_once(
        "corda2_step1", "corda2_step1_HC_dependencies", 4L,
        "supporting high-confidence directions with MC/NC dependencies"
      )
    } else if (identical(stage, "corda2_step2_1_MC_NC_dependencies")) {
      .rc_layer2_algorithm_once(
        "corda2_step2_1", "corda2_step2_1_MC_NC_dependencies", 5L,
        "measuring NC dependencies of remaining MC directions"
      )
    } else if (identical(stage, "corda2_step3_HC_OT_dependencies")) {
      .rc_layer2_algorithm_once(
        "corda2_step3", "corda2_step3_HC_OT_dependencies", 7L,
        "adding only OT reactions required by retained HC flux"
      )
    }
  }
  previous <- .rc_layer2_progress_state$inside_dependency
  .rc_layer2_progress_state$inside_dependency <- TRUE
  on.exit({ .rc_layer2_progress_state$inside_dependency <- previous }, add = TRUE)
  do.call(
    .rc_layer2_progress_original$.rc_corda2_dependency_assessment,
    args
  )
}

.rc_layer2_progress_capture(".rc_corda2_maximize_target")
.rc_corda2_maximize_target <- function(...) {
  task <- .rc_layer2_progress_state$current_task
  if (!is.null(task) && identical(task$route, "corda2") &&
      !isTRUE(.rc_layer2_progress_state$inside_dependency)) {
    .rc_layer2_algorithm_once(
      "corda2_step2_2", "corda2_step2_2_MC_feasibility", 6L,
      "promoting frequent NC dependencies and testing MC feasibility"
    )
  }
  do.call(.rc_layer2_progress_original$.rc_corda2_maximize_target, list(...))
}

.rc_layer2_progress_capture(".rc_corda_build_three_stage")
.rc_corda_build_three_stage <- function(...) {
  answer <- do.call(
    .rc_layer2_progress_original$.rc_corda_build_three_stage,
    list(...)
  )
  task <- .rc_layer2_progress_state$current_task
  if (!is.null(task) && identical(task$route, "corda2")) {
    required <- list(
      c("corda2_step1", "corda2_step1_HC_dependencies", 4L),
      c("corda2_step2_1", "corda2_step2_1_MC_NC_dependencies", 5L),
      c("corda2_step2_2", "corda2_step2_2_MC_feasibility", 6L),
      c("corda2_step3", "corda2_step3_HC_OT_dependencies", 7L)
    )
    for (item in required) {
      if (!exists(item[[1L]], envir = .rc_layer2_progress_state$algorithm_flags,
                  inherits = FALSE)) {
        .rc_layer2_algorithm_once(
          item[[1L]], item[[2L]], as.integer(item[[3L]]),
          "skipped because this confidence class had no candidate directions"
        )
      }
    }
  }
  answer
}

.rc_layer2_progress_capture(".rc_corda_core_closure")
.rc_corda_core_closure <- function(...) {
  task <- .rc_layer2_progress_state$current_task
  if (!is.null(task) && identical(task$route, "corda2")) {
    .rc_layer2_current_task_event(
      "corda2_target_closure", 8L,
      "validating retained core directions in the reconstructed GEM"
    )
  }
  do.call(.rc_layer2_progress_original$.rc_corda_core_closure, list(...))
}

# Full-GEM model construction has no cell-type partition.
.rc_layer2_progress_capture("rc_build_full_gem")
rc_build_full_gem <- function(...) {
  if (!isTRUE(.rc_layer2_progress_state$active) &&
      is.null(.rc_layer2_progress_state$current_task)) {
    return(do.call(.rc_layer2_progress_original$rc_build_full_gem, list(...)))
  }
  if (!is.null(.rc_layer2_progress_state$current_task)) {
    return(do.call(.rc_layer2_progress_original$rc_build_full_gem, list(...)))
  }
  args <- list(...)
  context <- .rc_layer2_task_context(
    "ALL", .rc_layer2_medium_id(args$medium_table), "full_gem"
  )
  parts_dir <- .rc_layer2_progress_dir_from_cache(
    .rc_layer2_progress_cache_dir_from_frames()
  )
  .rc_layer2_task_event(
    context, "full_gem_medium_bounds", 1L, 2L,
    "applying medium bounds without reaction deletion",
    scope = "structural", run_kind = "primary", parts_dir = parts_dir
  )
  answer <- do.call(.rc_layer2_progress_original$rc_build_full_gem, args)
  .rc_layer2_task_event(
    context, "full_gem_model_ready", 2L, 2L,
    detail = paste0("reactions=", ncol(answer$S)),
    scope = "structural", run_kind = "primary", status = "complete",
    parts_dir = parts_dir
  )
  answer
}

# Per-model vmax and Step 2 scoring phases.
.rc_layer2_progress_capture(".rc_build_microcompass_vmax_cache")
.rc_build_microcompass_vmax_cache <- function(...) {
  args <- list(...)
  contexts <- .rc_layer2_model_contexts(
    args$model_cache,
    mode = as.character(args$mode %||% "meta_module_gem")
  )
  parts_dir <- .rc_layer2_cache_progress_dir(args$model_cache)
  run_kind <- .rc_layer2_progress_state$run_kind %||% "primary"
  for (item in contexts) {
    .rc_layer2_task_event(
      item$context, "directional_vmax_start", 1L, 4L,
      detail = paste0("directional_targets=", item$n_targets),
      scope = "scoring", run_kind = run_kind, parts_dir = parts_dir
    )
  }
  answer <- do.call(
    .rc_layer2_progress_original$.rc_build_microcompass_vmax_cache,
    args
  )
  for (item in contexts) {
    .rc_layer2_task_event(
      item$context, "directional_vmax_complete", 2L, 4L,
      detail = paste0("directional_targets=", item$n_targets),
      scope = "scoring", run_kind = run_kind,
      status = "complete", parts_dir = parts_dir
    )
    .rc_layer2_task_event(
      item$context, "penalty_step2_start", 3L, 4L,
      "scoring matching metacells with the cached directional vmax",
      scope = "scoring", run_kind = run_kind, parts_dir = parts_dir
    )
  }
  if (identical(run_kind, "primary")) {
    .rc_layer2_overall_event(
      "primary_vmax_complete", 4L,
      detail = paste0(
        "directional target batches completed; tasks=",
        attr(answer, "parallel_tasks") %||% length(answer)
      )
    )
  }
  answer
}

.rc_layer2_progress_capture(".rc_run_celltype_microcompass_engine")
.rc_run_celltype_microcompass_engine <- function(...) {
  args <- list(...)
  answer <- do.call(
    .rc_layer2_progress_original$.rc_run_celltype_microcompass_engine,
    args
  )
  contexts <- .rc_layer2_model_contexts(
    answer$shared_model_cache, mode = "meta_module_gem"
  )
  parts_dir <- .rc_layer2_cache_progress_dir(answer$shared_model_cache)
  run_kind <- .rc_layer2_progress_state$run_kind %||% "primary"
  celltype_col <- as.character(args$celltype_col %||% "cell_type")
  unit_celltypes <- if (celltype_col %in% colnames(answer$unit_meta)) {
    as.character(answer$unit_meta[[celltype_col]])
  } else {
    character()
  }
  for (item in contexts) {
    .rc_layer2_task_event(
      item$context, "penalty_step2_complete", 4L, 4L,
      detail = paste0(
        "matching_units=", sum(unit_celltypes == item$context$cell_type)
      ),
      scope = "scoring", run_kind = run_kind,
      status = "complete", parts_dir = parts_dir
    )
  }
  answer
}

.rc_layer2_progress_capture(".rc_run_shared_full_gem_engine")
.rc_run_shared_full_gem_engine <- function(...) {
  answer <- do.call(
    .rc_layer2_progress_original$.rc_run_shared_full_gem_engine,
    list(...)
  )
  contexts <- .rc_layer2_model_contexts(
    answer$shared_model_cache, mode = "full_gem"
  )
  parts_dir <- .rc_layer2_cache_progress_dir(answer$shared_model_cache)
  run_kind <- .rc_layer2_progress_state$run_kind %||% "primary"
  for (item in contexts) {
    .rc_layer2_task_event(
      item$context, "penalty_step2_complete", 4L, 4L,
      detail = paste0("units=", nrow(answer$unit_meta)),
      scope = "scoring", run_kind = run_kind,
      status = "complete", parts_dir = parts_dir
    )
  }
  answer
}

# Structural cache completion advances the controller-owned overall bar.
.rc_layer2_progress_capture(".rc_build_celltype_medium_union_gem_cache")
.rc_build_celltype_medium_union_gem_cache <- function(...) {
  answer <- do.call(
    .rc_layer2_progress_original$.rc_build_celltype_medium_union_gem_cache,
    list(...)
  )
  if (identical(.rc_layer2_progress_state$run_kind, "primary")) {
    summary <- attr(answer, "summary")
    .rc_layer2_overall_event(
      "structural_models_complete", 3L,
      detail = paste0(
        "cell_type_x_medium_models=",
        if (is.data.frame(summary)) nrow(summary) else length(unique(vapply(
          answer, function(entry) as.character(entry$file), character(1)
        )))
      )
    )
  }
  answer
}

.rc_layer2_progress_capture("rc_build_full_gem_cache")
rc_build_full_gem_cache <- function(...) {
  answer <- do.call(
    .rc_layer2_progress_original$rc_build_full_gem_cache,
    list(...)
  )
  if (identical(.rc_layer2_progress_state$run_kind, "primary")) {
    summary <- attr(answer, "summary")
    .rc_layer2_overall_event(
      "structural_models_complete", 3L,
      detail = paste0(
        "full_gem_medium_models=",
        if (is.data.frame(summary)) nrow(summary) else 1L
      )
    )
  }
  answer
}

# Final comparison and validation phases.
.rc_layer2_progress_capture(".rc_layer2_comparison_table")
.rc_layer2_comparison_table <- function(...) {
  answer <- do.call(
    .rc_layer2_progress_original$.rc_layer2_comparison_table,
    list(...)
  )
  .rc_layer2_overall_event(
    "comparison_table_complete", 9L,
    detail = paste0("comparison_rows=", nrow(answer))
  )
  answer
}

.rc_layer2_progress_capture(".rc_validate_layer2_stage")
.rc_validate_layer2_stage <- function(...) {
  answer <- do.call(
    .rc_layer2_progress_original$.rc_validate_layer2_stage,
    list(...)
  )
  .rc_layer2_overall_event(
    "layer2_validation_complete", 10L,
    "Layer 2 schemas, ordering and shared-model contracts validated"
  )
  answer
}
