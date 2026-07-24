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
  .rc_progress_update(state, 0L, "started")
  state
}

.rc_progress_update <- function(state, current, detail = NULL) {
  if (is.null(state) || !is.environment(state) || !isTRUE(state$enabled)) {
    return(invisible(state))
  }
  current <- max(0L, min(state$total, as.integer(current[[1L]])))
  state$current <- current
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
  message(sprintf(
    "%s [%s] %d/%d%s",
    state$label, bar, current, state$total, suffix
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

.rc_step_monitor_start <- function(
    stage, outdir, progress = TRUE, total_parts = 1L) {
  progress <- .rc_progress_enabled(progress)
  monitor <- new.env(parent = emptyenv())
  monitor$timer <- .rc_timing_start(stage)
  monitor$outdir <- outdir
  monitor$old_progress_option <- options(RegCompassR.progress = progress)
  monitor$progress <- .rc_progress_new(
    total_parts, paste0("RegCompass ", stage), progress
  )
  monitor$finished <- FALSE
  monitor$option_restored <- FALSE
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
    monitor$timer, status = status, outdir = monitor$outdir, details = details
  )
  monitor$finished <- TRUE
  .rc_progress_done(monitor$progress, status)
  .rc_restore_monitor_progress_option(monitor)
  if (is.list(value)) value$timing <- timing
  value
}

.rc_step_monitor_fail <- function(monitor) {
  if (!is.null(monitor) && is.environment(monitor) &&
      !isTRUE(monitor$finished)) {
    .rc_timing_finish(
      monitor$timer, status = "error", outdir = monitor$outdir
    )
    .rc_progress_done(monitor$progress, "error")
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
