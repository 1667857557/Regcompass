`%||%` <- function(x, y) if (is.null(x)) y else x

.rc_layer2_completion_context <- new.env(parent = emptyenv())
.rc_layer2_completion_context$active <- FALSE
.rc_layer2_completion_context$model_completion <- "fastcore"
.rc_layer2_completion_context$corda_options <- NULL
.rc_is_corda2_options <- function(x) FALSE

source("R/layer2_parallel_runtime.R")
source("R/celltype_union_gem_cache.R")

.rc_validate_shared_medium <- function(x) x
.rc_bind_frames_fill <- function(frames) {
  frames <- frames[vapply(frames, is.data.frame, logical(1))]
  if (!length(frames)) return(data.frame())
  columns <- unique(unlist(lapply(frames, colnames), use.names = FALSE))
  frames <- lapply(frames, function(frame) {
    missing <- setdiff(columns, colnames(frame))
    for (name in missing) frame[[name]] <- NA
    frame[columns]
  })
  answer <- do.call(rbind, frames)
  rownames(answer) <- NULL
  answer
}

observed_nested <- logical()
.rc_complete_celltype_medium_union_gem <- function(
    gem, reaction_membership, core_reactions, cell_type, medium_table = NULL,
    target_direction = "both", solver = "highs", time_limit = 300,
    fastcore_epsilon = 1e-4, max_support_reactions = 2000,
    strict = TRUE) {
  observed_nested <<- c(
    observed_nested,
    isTRUE(.rc_layer2_parallel_context$nested_serial)
  )
  stopifnot(
    is.data.frame(reaction_membership),
    is.data.frame(core_reactions),
    nrow(core_reactions) == 1L
  )
  list(
    S = matrix(
      c(1, -1), nrow = 1L,
      dimnames = list("M", c("Rcore", "Rsupport"))
    ),
    lb = c(Rcore = 0, Rsupport = 0),
    ub = c(Rcore = 1000, Rsupport = 1000),
    target_directions = data.frame(
      reaction_id = "Rcore",
      target_direction = "forward",
      stringsAsFactors = FALSE
    ),
    build_params = list(
      n_celltype_biological_reactions = nrow(reaction_membership),
      n_celltype_fastcore_support_reactions = 1L
    ),
    target_status = "ok"
  )
}

parallel_calls <- new.env(parent = emptyenv())
parallel_calls$n <- 0L
parallel_calls$tasks <- 0L
rc_parallel_lapply <- function(X, FUN, BPPARAM = NULL, ...) {
  parallel_calls$n <- parallel_calls$n + 1L
  parallel_calls$tasks <- length(X)
  lapply(X, FUN)
}

# Exercise the outer-parallel branch without requiring BiocParallel in this
# lightweight source regression. Package CI separately parses the real backend.
.rc_layer2_tune_task_bpparam <- function(BPPARAM, n_tasks) BPPARAM
.rc_layer2_pool_workers <- function(BPPARAM) {
  if (identical(BPPARAM, FALSE)) 1L else 4L
}
.rc_layer2_should_outer_parallel <- function(n_tasks, pool_workers) {
  n_tasks > 1L && pool_workers > 1L
}

membership <- data.frame(
  cell_type = rep(c("A", "B"), each = 2L),
  reaction_id = rep(c("Rcore", "Rsupport"), times = 2L),
  stringsAsFactors = FALSE
)
core <- data.frame(
  cell_type = c("A", "B"),
  reaction_id = "Rcore",
  is_core = TRUE,
  stringsAsFactors = FALSE
)
media <- data.frame(
  medium_scenario_id = c("m1", "m2"),
  .no_constraints = TRUE,
  stringsAsFactors = FALSE
)

outdir <- tempfile("fastcore-parallel-cache-")
dir.create(outdir, recursive = TRUE)
.rc_layer2_parallel_context$active <- TRUE
.rc_layer2_parallel_context$parallel <- TRUE
.rc_layer2_parallel_context$BPPARAM <- structure(list(), class = "fake_param")
.rc_layer2_parallel_context$nested_serial <- FALSE

cache <- .rc_build_celltype_medium_union_gem_cache(
  gem = list(dummy = TRUE),
  reaction_membership = membership,
  core_reactions = core,
  medium_scenarios = media,
  celltype_col = "cell_type",
  cache_dir = outdir,
  strict = TRUE
)

summary <- attr(cache, "summary")
stopifnot(
  identical(parallel_calls$n, 1L),
  identical(parallel_calls$tasks, 4L),
  length(cache) == 4L,
  is.data.frame(summary),
  nrow(summary) == 4L,
  identical(attr(cache, "completion_method"), "fastcore"),
  identical(attr(cache, "structural_parallel_workers"), 4L),
  identical(attr(cache, "structural_parallel_tasks"), 4L),
  isTRUE(attr(cache, "structural_dynamic_task_scheduling")),
  grepl("outer_parallel_fastcore", attr(cache, "fastcore_parallel_task")),
  grepl("drop task-local", attr(cache, "fastcore_worker_cleanup")),
  all(observed_nested),
  !isTRUE(.rc_layer2_parallel_context$nested_serial),
  all(file.exists(summary$file)),
  !anyDuplicated(names(cache)),
  !length(list.files(outdir, pattern = "[.]tmp_", recursive = TRUE))
)

for (file in summary$file) {
  model <- readRDS(file)
  stopifnot(
    isTRUE(model$is_union_gem),
    identical(
      model$union_gem_scope,
      "one_cell_type_one_medium_shared_across_conditions_and_matching_metacells"
    )
  )
}

# Serial mode must bypass the parallel dispatcher while retaining the same
# cache/result contract.
parallel_calls$n <- 0L
parallel_calls$tasks <- 0L
observed_nested <- logical()
serial_dir <- tempfile("fastcore-serial-cache-")
dir.create(serial_dir, recursive = TRUE)
.rc_layer2_parallel_context$parallel <- FALSE
.rc_layer2_parallel_context$BPPARAM <- FALSE
serial_cache <- .rc_build_celltype_medium_union_gem_cache(
  gem = list(dummy = TRUE),
  reaction_membership = membership,
  core_reactions = core,
  medium_scenarios = media[1L, , drop = FALSE],
  celltype_col = "cell_type",
  cache_dir = serial_dir,
  strict = TRUE
)
stopifnot(
  identical(parallel_calls$n, 0L),
  length(serial_cache) == 2L,
  identical(attr(serial_cache, "structural_parallel_workers"), 1L),
  identical(attr(serial_cache, "structural_parallel_tasks"), 2L),
  !isTRUE(attr(serial_cache, "structural_dynamic_task_scheduling")),
  all(!observed_nested)
)

unlink(outdir, recursive = TRUE, force = TRUE)
unlink(serial_dir, recursive = TRUE, force = TRUE)
message("FASTCORE parallel cache and worker-cleanup contract passed.")
