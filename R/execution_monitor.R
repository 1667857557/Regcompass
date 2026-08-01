.rc_progress_enabled <- function(
    progress = getOption("RegCompassR.progress", TRUE)) {
  if (is.null(progress)) progress <- TRUE
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("`progress` must be TRUE or FALSE.", call. = FALSE)
  }
  isTRUE(progress)
}

.rc_progress_new <- function(total, label, progress = TRUE) {
  total <- suppressWarnings(as.integer(total[[1L]]))
  if (!is.finite(total) || total < 1L) total <- 1L
  state <- new.env(parent = emptyenv())
  state$total <- total
  state$current <- 0L
  state$label <- as.character(label[[1L]])
  state$enabled <- .rc_progress_enabled(progress)
  state$started_elapsed <- unname(proc.time()[["elapsed"]])
  .rc_progress_update(state, 0L, "started")
  state
}

.rc_progress_update <- function(state, current, detail = NULL) {
  if (is.null(state) || !is.environment(state)) {
    return(invisible(state))
  }
  current <- max(0L, min(state$total, as.integer(current[[1L]])))
  state$current <- current
  if (!isTRUE(state$enabled)) return(invisible(state))
  width <- 24L
  filled <- floor(width * current / state$total)
  cursor <- as.integer(filled < width)
  bar <- paste0(
    strrep("=", filled),
    if (cursor) ">" else "",
    strrep(" ", max(0L, width - filled - cursor))
  )
  suffix <- if (is.null(detail) || !nzchar(as.character(detail[[1L]]))) {
    ""
  } else {
    paste0(" ", as.character(detail[[1L]]))
  }
  elapsed <- max(
    0,
    unname(proc.time()[["elapsed"]]) - as.numeric(state$started_elapsed)
  )
  percent <- 100 * current / state$total
  message(sprintf(
    "%s [%s] %d/%d (%5.1f%%) elapsed=%s%s",
    state$label, bar, current, state$total, percent,
    .rc_format_elapsed(elapsed), suffix
  ))
  invisible(state)
}

.rc_progress_done <- function(state, detail = "complete") {
  if (!is.null(state) && is.environment(state)) {
    .rc_progress_update(state, state$total, detail)
  }
  invisible(state)
}

.rc_format_elapsed <- function(seconds) {
  seconds <- max(0, as.numeric(seconds[[1L]]))
  hours <- floor(seconds / 3600)
  minutes <- floor((seconds %% 3600) / 60)
  secs <- seconds %% 60
  sprintf("%02d:%02d:%06.3f", hours, minutes, secs)
}

.rc_progress_context_string <- function(context) {
  if (is.null(context) || !length(context)) return("")
  if (!is.list(context) || is.null(names(context)) || any(!nzchar(names(context)))) {
    stop("Progress event context must be a named list.", call. = FALSE)
  }
  values <- vapply(names(context), function(name) {
    value <- context[[name]]
    if (is.null(value) || !length(value)) value <- NA_character_
    value <- paste(as.character(value), collapse = ",")
    value <- gsub("[\t\r\n;]+", " ", value)
    paste0(name, "=", value)
  }, character(1))
  paste(values, collapse = ";")
}

.rc_step_monitor_event <- function(
    monitor, phase, detail = NULL, current = NULL, total = NULL,
    context = list(), status = "running", emit = TRUE) {
  if (is.null(monitor) || !is.environment(monitor)) {
    return(invisible(monitor))
  }
  phase <- as.character(phase[[1L]])
  status <- as.character(status[[1L]])
  if (!nzchar(phase) || !nzchar(status)) {
    stop("Progress event phase and status must be non-empty.", call. = FALSE)
  }
  if (is.null(current)) current <- monitor$progress$current
  if (is.null(total)) total <- monitor$progress$total
  current <- max(0L, as.integer(current[[1L]]))
  total <- max(1L, as.integer(total[[1L]]))
  current <- min(current, total)
  elapsed <- max(
    0,
    unname(proc.time()[["elapsed"]]) -
      as.numeric(monitor$timer$elapsed_start)
  )
  context_text <- .rc_progress_context_string(context)
  detail_text <- if (is.null(detail) || !length(detail)) "" else {
    gsub("[\t\r\n]+", " ", as.character(detail[[1L]]))
  }
  monitor$event_sequence <- as.integer(monitor$event_sequence %||% 0L) + 1L
  monitor$last_phase <- phase
  monitor$last_detail <- detail_text
  row <- data.frame(
    sequence = monitor$event_sequence,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3%z"),
    stage = as.character(monitor$timer$stage),
    phase = phase,
    status = status,
    current = current,
    total = total,
    percent = 100 * current / total,
    elapsed_seconds = elapsed,
    elapsed_hms = .rc_format_elapsed(elapsed),
    detail = detail_text,
    context = context_text,
    stringsAsFactors = FALSE
  )
  if (!is.null(monitor$progress_log)) {
    append <- file.exists(monitor$progress_log)
    utils::write.table(
      row,
      file = monitor$progress_log,
      sep = "\t", quote = FALSE, row.names = FALSE,
      col.names = !append, append = append
    )
  }
  if (isTRUE(emit)) {
    console_detail <- paste0(
      "phase=", phase,
      if (nzchar(detail_text)) paste0(" | ", detail_text) else "",
      if (nzchar(context_text)) paste0(" | ", context_text) else ""
    )
    .rc_progress_update(monitor$progress, current, console_detail)
  } else {
    monitor$progress$current <- current
  }
  invisible(row)
}

.rc_timing_start <- function(stage) {
  list(
    stage = as.character(stage[[1L]]),
    started_at = Sys.time(),
    elapsed_start = unname(proc.time()[["elapsed"]])
  )
}

.rc_timing_finish <- function(
    timer, status = "success", outdir = NULL, details = NULL) {
  finished_at <- Sys.time()
  elapsed_seconds <- max(
    0,
    unname(proc.time()[["elapsed"]]) - as.numeric(timer$elapsed_start)
  )
  row <- data.frame(
    stage = as.character(timer$stage),
    status = as.character(status),
    started_at = format(timer$started_at, "%Y-%m-%dT%H:%M:%S%z"),
    finished_at = format(finished_at, "%Y-%m-%dT%H:%M:%S%z"),
    elapsed_seconds = elapsed_seconds,
    elapsed_hms = .rc_format_elapsed(elapsed_seconds),
    os_type = .Platform$OS.type,
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    stringsAsFactors = FALSE
  )
  if (!is.null(details) && is.list(details) && length(details)) {
    for (name in names(details)) {
      value <- details[[name]]
      if (length(value) != 1L) value <- paste(value, collapse = ";")
      row[[name]] <- value
    }
  }
  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    utils::write.table(
      row,
      file = file.path(outdir, "step_timing.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
    )
  }
  row
}

.rc_step_final_artifact <- function(stage, outdir) {
  filename <- switch(
    as.character(stage[[1L]]),
    grn = "step_grn.rds",
    metacells = "step_metacells.rds",
    meta_modules = "step_meta_modules.rds",
    layer1 = "step_layer1.rds",
    layer2 = "step_layer2.rds",
    results = "regcompass_result.rds",
    target_union = "step_target_union.rds",
    NULL
  )
  if (is.null(filename)) return(NULL)
  file.path(outdir, filename)
}

.rc_step_artifact_signature <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    return(list(exists = FALSE, size = NA_real_, mtime = NA_real_, ctime = NA_real_))
  }
  info <- file.info(path)
  list(
    exists = TRUE,
    size = as.numeric(info$size[[1L]]),
    mtime = as.numeric(info$mtime[[1L]]),
    ctime = as.numeric(info$ctime[[1L]])
  )
}

.rc_step_artifact_committed <- function(path, before) {
  after <- .rc_step_artifact_signature(path)
  if (!isTRUE(after$exists)) return(FALSE)
  if (!is.list(before) || !isTRUE(before$exists)) return(TRUE)
  !identical(after$size, before$size) ||
    (!is.na(after$mtime) && !is.na(before$mtime) && after$mtime > before$mtime) ||
    (!is.na(after$ctime) && !is.na(before$ctime) && after$ctime > before$ctime)
}

.rc_step_monitor_start <- function(
    stage, outdir, progress = TRUE, total_parts = 1L) {
  progress <- .rc_progress_enabled(progress)
  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }
  monitor <- new.env(parent = emptyenv())
  monitor$timer <- .rc_timing_start(stage)
  monitor$outdir <- outdir
  monitor$old_progress_option <- options(RegCompassR.progress = progress)
  monitor$progress <- .rc_progress_new(
    total_parts, paste0("RegCompass ", stage), progress
  )
  monitor$progress_log <- if (is.null(outdir)) NULL else {
    file.path(outdir, "step_progress.tsv")
  }
  if (!is.null(monitor$progress_log) && file.exists(monitor$progress_log)) {
    unlink(monitor$progress_log)
  }
  monitor$event_sequence <- 0L
  monitor$final_artifact <- .rc_step_final_artifact(stage, outdir)
  monitor$artifact_before <- .rc_step_artifact_signature(
    monitor$final_artifact
  )
  monitor$finish_requested <- FALSE
  monitor$finish_status <- "success"
  monitor$finish_details <- NULL
  monitor$finished <- FALSE
  monitor$option_restored <- FALSE
  .rc_step_monitor_event(
    monitor,
    phase = "stage_start",
    detail = "stage monitor initialized",
    current = 0L,
    status = "running",
    context = list(progress_enabled = progress),
    emit = FALSE
  )
  monitor
}

.rc_restore_monitor_progress_option <- function(monitor) {
  if (!is.null(monitor) && is.environment(monitor) &&
      !isTRUE(monitor$option_restored)) {
    do.call(options, monitor$old_progress_option)
    monitor$option_restored <- TRUE
  }
  invisible(NULL)
}

.rc_step_monitor_finish <- function(
    value, monitor, status = "success", details = NULL) {
  if (is.null(monitor) || !is.environment(monitor)) return(value)
  timing <- .rc_timing_finish(
    monitor$timer, status = status, outdir = NULL, details = details
  )
  finalizing_current <- max(0L, monitor$progress$total - 1L)
  .rc_step_monitor_event(
    monitor,
    phase = "final_artifact",
    detail = "writing final artifacts",
    current = finalizing_current,
    status = "running"
  )
  if (is.null(monitor$final_artifact)) {
    .rc_timing_finish(
      monitor$timer, status = status, outdir = monitor$outdir,
      details = details
    )
    .rc_step_monitor_event(
      monitor, phase = "stage_complete", detail = status,
      current = monitor$progress$total, status = status, emit = FALSE
    )
    monitor$finished <- TRUE
    .rc_progress_done(monitor$progress, status)
    .rc_restore_monitor_progress_option(monitor)
  } else {
    monitor$finish_requested <- TRUE
    monitor$finish_status <- status
    monitor$finish_details <- details
  }
  if (is.list(value)) value$timing <- timing
  value
}

.rc_step_monitor_fail <- function(monitor) {
  if (!is.null(monitor) && is.environment(monitor) &&
      !isTRUE(monitor$finished)) {
    artifact_committed <- isTRUE(monitor$finish_requested) &&
      .rc_step_artifact_committed(
        monitor$final_artifact, monitor$artifact_before
      )
    status <- if (artifact_committed) monitor$finish_status else "error"
    details <- if (artifact_committed) monitor$finish_details else NULL
    .rc_timing_finish(
      monitor$timer, status = status, outdir = monitor$outdir,
      details = details
    )
    diagnostic <- paste0(
      "stage stopped before completion; last_phase=",
      monitor$last_phase %||% "stage_start",
      if (nzchar(monitor$last_detail %||% "")) {
        paste0("; last_detail=", monitor$last_detail)
      } else "",
      ". Inspect step_progress.tsv and step_timing.tsv for diagnosis."
    )
    .rc_step_monitor_event(
      monitor,
      phase = if (artifact_committed) "stage_complete" else "stage_error",
      detail = if (artifact_committed) status else diagnostic,
      current = if (artifact_committed) {
        monitor$progress$total
      } else {
        monitor$progress$current
      },
      status = status,
      emit = !artifact_committed
    )
    .rc_progress_done(monitor$progress, status)
    monitor$finished <- TRUE
  }
  .rc_restore_monitor_progress_option(monitor)
  invisible(NULL)
}

.rc_write_execution_timing <- function(timing, outdir) {
  if (!is.data.frame(timing)) {
    stop("`timing` must be a data frame.", call. = FALSE)
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    timing,
    file = file.path(outdir, "00_execution_timing.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
  )
  invisible(timing)
}
