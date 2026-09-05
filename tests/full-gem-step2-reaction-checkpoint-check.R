suppressPackageStartupMessages({
  library(Matrix)
  library(highs)
})

`%||%` <- function(x, y) if (is.null(x)) y else x
.rc_as_dgCMatrix <- function(x) methods::as(x, "dgCMatrix")

rc_align_bound <- function(x, reactions, default, name) {
  if (is.null(x)) {
    value <- rep(default, length(reactions))
    names(value) <- reactions
    return(value)
  }
  if (!is.null(names(x))) x <- x[reactions]
  x <- as.numeric(x)
  if (length(x) != length(reactions) || anyNA(x)) {
    stop("invalid ", name, " bounds", call. = FALSE)
  }
  names(x) <- reactions
  x
}

source("R/00_utils.R")
source("R/lp_solver.R")
source("R/microcompass_vmax_cache.R")
source("R/layer2_step2_directional_prepare.R")
source("R/celltype_microcompass_reaction_parallel.R")
source("R/layer2_step2_model_batch.R")

.rc_atomic_save_rds <- function(x, file) {
  saveRDS(x, file)
  invisible(file)
}
.rc_bind_frames_fill <- function(frames) {
  frames <- frames[vapply(frames, is.data.frame, logical(1))]
  if (!length(frames)) return(data.frame())
  do.call(rbind, frames)
}
.rc_microcompass_object_checksum <- function(x) {
  raw <- serialize(x, NULL, version = 2L)
  paste0(length(raw), "-", sum(as.integer(raw)))
}

root <- tempfile("full-gem-step2-model-batch-ci-")
dir.create(root, recursive = TRUE)
on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

# Full-GEM forward-only model deliberately exercises the zero-length reverse
# direction edge case that previously produced a spurious "::reverse" ID.
S <- Matrix::Matrix(
  matrix(
    c(1, -1, -1), nrow = 1,
    dimnames = list("M", c("UP", "T1", "T2"))
  ),
  sparse = TRUE
)
lb <- c(UP = 0, T1 = 0, T2 = 0)
ub <- c(UP = 10, T1 = 10, T2 = 10)
reactions <- c("UP", "T1", "T2")
units <- c("mc1", "mc2", "mc3")
row_ids <- c(
  "reaction=T1::direction=forward::medium=toy::condition=all",
  "reaction=T2::direction=forward::medium=toy::condition=all"
)
entries <- stats::setNames(lapply(c("T1", "T2"), function(target) {
  list(
    reaction_id = target,
    target_direction = "forward",
    medium_scenario = "toy",
    condition = "all"
  )
}), row_ids)

penalty_matrix <- matrix(
  c(
    0.20, 0.50, 0.90,
    0.25, 0.80, 0.40,
    0.70, 0.10, 0.60
  ),
  nrow = 3,
  dimnames = list(reactions, units)
)
control_penalty <- penalty_matrix
control_penalty[, "mc1"] <- penalty_matrix[, "mc1"]
control_penalty[, "mc2"] <- c(0.45, 0.35, 0.75)
control_penalty[, "mc3"] <- c(0.10, 0.95, 0.25)
control_identical <- c(mc1 = TRUE, mc2 = FALSE, mc3 = FALSE)

vmax_values <- stats::setNames(lapply(c("T1", "T2"), function(target) {
  .rc_step2_compact_vmax_value(rc_compass_vmax_directional(
    S, lb, ub, target, direction = "forward", solver = "highs"
  ))
}), row_ids)
stopifnot(all(vapply(vmax_values, function(x) isTRUE(x$feasible), logical(1))))

payload <- list(
  schema_version = "regcompass_full_gem_step2_compact_payload_v1",
  model = list(
    S = S,
    lb = lb,
    ub = ub,
    target_status = "not_prechecked",
    closure_diagnostics = data.frame(),
    file_checksum = "full-gem-shared-toy"
  ),
  reactions = reactions,
  units = units,
  penalty = penalty_matrix,
  control_penalty = control_penalty,
  control_identical = control_identical,
  entries = entries,
  vmax = vmax_values,
  omega = 0.95,
  solver = "highs",
  flux_threshold = 1e-8
)

payload_file <- file.path(root, "payload.rds")
checkpoint_dir <- file.path(root, "checkpoints")
dir.create(checkpoint_dir)
saveRDS(payload, payload_file)

checkpoints <- .rc_full_gem_step2_reaction_batch_worker(list(
  payload_file = payload_file,
  row_ids = row_ids,
  checkpoint_dir = checkpoint_dir
))
stopifnot(
  length(checkpoints) == length(row_ids),
  all(file.exists(checkpoints))
)

observed_rows <- character()
for (checkpoint in checkpoints) {
  observed <- readRDS(checkpoint)
  row_id <- observed$row_id
  observed_rows <- c(observed_rows, row_id)
  entry <- entries[[row_id]]

  stopifnot(
    identical(observed$units, units),
    observed$engine_metrics$n_solves == length(units),
    isTRUE(observed$engine_metrics$shared_model_batch_engine),
    observed$engine_metrics$batch_objective_change_events <= length(units),
    observed$engine_metrics$batch_target_switches ==
      length(row_ids) * length(units),
    observed$control$engine_metrics$n_solves == sum(!control_identical),
    all(observed$diagnostics$step2_solver_reused_across_targets),
    all(observed$diagnostics$step2_traversal == "unit_then_directional_target"),
    identical(
      observed$control$diagnostics$objective_identical_to_primary,
      unname(control_identical)
    )
  )

  for (one_unit in units) {
    oracle_primary <- rc_compass_two_step_lp_directional(
      S, lb, ub, entry$reaction_id, penalty_matrix[, one_unit],
      target_direction = entry$target_direction,
      omega = 0.95, solver = "highs"
    )
    oracle_control <- rc_compass_two_step_lp_directional(
      S, lb, ub, entry$reaction_id, control_penalty[, one_unit],
      target_direction = entry$target_direction,
      omega = 0.95, solver = "highs"
    )
    stopifnot(
      identical(observed$feasible[[one_unit]], oracle_primary$feasible),
      isTRUE(all.equal(
        observed$vmax[[one_unit]], oracle_primary$vmax,
        tolerance = 1e-10
      )),
      isTRUE(all.equal(
        observed$penalty[[one_unit]], oracle_primary$penalty,
        tolerance = 1e-9
      )),
      identical(
        observed$control$feasible[[one_unit]], oracle_control$feasible
      ),
      isTRUE(all.equal(
        observed$control$penalty[[one_unit]], oracle_control$penalty,
        tolerance = 1e-9
      ))
    )
  }
}
stopifnot(setequal(observed_rows, row_ids))

# Preserve target/reaction outer parallelism. One shared full GEM with enough
# targets must still expose the requested 80 independent target batches.
parallel_rows <- paste0("directional_row_", seq_len(200L))
parallel_model_keys <- stats::setNames(
  rep("shared_full_gem", length(parallel_rows)), parallel_rows
)
parallel_batches <- .rc_step2_model_batches(
  parallel_model_keys, workers = 80L
)
parallel_assigned <- unlist(
  lapply(parallel_batches, `[[`, "row_ids"), use.names = FALSE
)
stopifnot(
  length(parallel_batches) == 80L,
  length(parallel_assigned) == length(parallel_rows),
  !anyDuplicated(parallel_assigned),
  setequal(parallel_assigned, parallel_rows)
)

cat("Full-GEM exact model-batch Step 2 numerical regression passed.\n")
