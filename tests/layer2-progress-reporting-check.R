`%||%` <- function(x, y) if (is.null(x)) y else x
.rc_bind_frames_fill <- function(frames) {
  frames <- frames[vapply(frames, is.data.frame, logical(1))]
  if (!length(frames)) return(data.frame())
  columns <- unique(unlist(lapply(frames, colnames), use.names = FALSE))
  aligned <- lapply(frames, function(frame) {
    missing <- setdiff(columns, colnames(frame))
    for (name in missing) frame[[name]] <- NA
    frame[columns]
  })
  do.call(rbind, aligned)
}

source("R/execution_monitor.R")
source("R/layer2_corda_pool_lifecycle.R")

probe_context <- .rc_layer2_task_context(
  cell_type = "Astrocyte",
  medium_scenario = "normal_human_plasma",
  route = "corda2"
)
probe_previous <- .rc_layer2_task_push(
  probe_context, route = "corda2", total = 9L
)
probe_dir <- .rc_layer2_progress_state$current_task$parts_dir
stopifnot(
  is.character(probe_dir),
  length(probe_dir) == 1L,
  !is.na(probe_dir),
  nzchar(probe_dir),
  dir.exists(probe_dir)
)
.rc_layer2_task_pop(probe_previous)
unlink(probe_dir, recursive = TRUE, force = TRUE)
.rc_layer2_progress_reset()

outdir <- tempfile("layer2-progress-")
dir.create(outdir, recursive = TRUE)
monitor <- .rc_step_monitor_start(
  "layer2", outdir = outdir, progress = FALSE
)
stopifnot(
  identical(monitor$progress$total, 12L),
  isTRUE(.rc_layer2_progress_state$active),
  identical(.rc_layer2_progress_state$owner_pid, Sys.getpid())
)

.rc_layer2_overall_event(
  "layer2_plan_ready", 1L,
  "two cell types x one medium",
  context = list(cell_types = 2L, media = 1L)
)
stopifnot(identical(monitor$progress$current, 1L))

task <- .rc_layer2_task_context(
  cell_type = "Astrocyte",
  medium_scenario = "normal_human_plasma",
  route = "corda2"
)
.rc_layer2_task_event(
  task, "task_start", 1L, 9L,
  "task inputs selected", parts_dir = .rc_layer2_progress_state$parts_dir,
  emit = FALSE
)
.rc_layer2_task_event(
  task, "corda2_step1_HC_dependencies", 4L, 9L,
  "supporting HC directions", parts_dir = .rc_layer2_progress_state$parts_dir,
  emit = FALSE
)
.rc_layer2_task_event(
  task, "corda2_model_ready", 9L, 9L,
  "included_reactions=42", status = "complete",
  parts_dir = .rc_layer2_progress_state$parts_dir, emit = FALSE
)

scoring <- .rc_layer2_task_context(
  cell_type = "Astrocyte",
  medium_scenario = "normal_human_plasma",
  route = "corda2"
)
.rc_layer2_task_event(
  scoring, "directional_vmax_complete", 2L, 4L,
  "directional_targets=8", scope = "scoring",
  parts_dir = .rc_layer2_progress_state$parts_dir, emit = FALSE
)
.rc_layer2_collect_task_progress(outdir)

progress_file <- file.path(outdir, "layer2_task_progress.tsv")
status_file <- file.path(outdir, "layer2_task_status.tsv")
stopifnot(file.exists(progress_file), file.exists(status_file))
events <- utils::read.delim(progress_file, stringsAsFactors = FALSE)
status <- utils::read.delim(status_file, stringsAsFactors = FALSE)
stopifnot(
  nrow(events) == 4L,
  all(c(
    "run_kind", "scope", "cell_type", "medium_scenario", "route",
    "phase", "step", "total", "percent", "elapsed_hms"
  ) %in% colnames(events)),
  any(events$phase == "corda2_step1_HC_dependencies"),
  any(events$scope == "scoring"),
  nrow(status) == 2L
)

cache <- list(
  a = list(
    file = file.path(outdir, "model_a.rds"),
    cell_type = "Astrocyte",
    medium_scenario = "plasma",
    build_strategy = "celltype_medium_original_matlab_corda2"
  ),
  b = list(
    file = file.path(outdir, "model_a.rds"),
    cell_type = "Astrocyte",
    medium_scenario = "plasma",
    build_strategy = "celltype_medium_original_matlab_corda2"
  ),
  c = list(
    file = file.path(outdir, "model_b.rds"),
    cell_type = "Microglia",
    medium_scenario = "plasma",
    build_strategy = "celltype_medium_union_gem"
  )
)
contexts <- .rc_layer2_model_contexts(cache, mode = "meta_module_gem")
stopifnot(
  length(contexts) == 2L,
  identical(contexts[[1L]]$context$route, "corda2"),
  identical(contexts[[1L]]$n_targets, 2L),
  identical(contexts[[2L]]$context$route, "fastcore")
)

.rc_layer2_progress_reset()
stopifnot(!isTRUE(.rc_layer2_progress_state$active))
unlink(outdir, recursive = TRUE, force = TRUE)
message("Layer 2 hierarchical progress reporting contract passed.")
