`%||%` <- function(x, y) if (is.null(x)) y else x

suppressPackageStartupMessages({
  library(Matrix)
  library(highs)
})

source("R/00_utils.R", local = FALSE)

rc_align_bound <- function(x, reactions, default, name) {
  if (is.null(x)) {
    value <- rep(default, length(reactions))
    names(value) <- reactions
    return(value)
  }
  if (!is.null(names(x))) {
    missing <- setdiff(reactions, names(x))
    if (length(missing)) stop(name, " is missing reactions: ", paste(missing, collapse = ", "))
    x <- x[reactions]
  }
  x <- as.numeric(x)
  if (length(x) != length(reactions) || anyNA(x)) stop(name, " does not align with reactions")
  names(x) <- reactions
  x
}

source("R/lp_solver.R", local = FALSE)
source("R/microcompass_vmax_cache.R", local = FALSE)
source("R/layer2_step2_directional_prepare.R", local = FALSE)
source("R/celltype_microcompass_reaction_parallel.R", local = FALSE)
source("R/layer2_step2_model_batch.R", local = FALSE)
source("R/layer2_step2_model_batch_streaming.R", local = FALSE)

.rc_atomic_save_rds <- function(x, file) {
  saveRDS(x, file)
  invisible(file)
}
.rc_microcompass_object_checksum <- function(x) {
  raw <- serialize(x, NULL, version = 2L)
  paste0(length(raw), "-", sum(as.integer(raw)))
}

# Engineering memory bound: retained target-by-unit result state cannot grow
# linearly with an arbitrarily large worker target batch.
stopifnot(
  .rc_step2_stream_target_chunk_size(1000L, 1000L) <= 16L,
  .rc_step2_stream_target_chunk_size(1000L, 1000L) * 1000L <= 100000L,
  .rc_step2_stream_target_chunk_size(200000L, 1000L) == 1L
)

# Truthful diagnostic contract.
stopifnot(
  .rc_step2_stream_cross_target_reuse(
    "highs_persistent_cpp", 2L, TRUE, 1L
  ),
  !.rc_step2_stream_cross_target_reuse(
    "one_shot", 2L, TRUE, 1L
  ),
  !.rc_step2_stream_cross_target_reuse(
    "highs_persistent_cpp", 1L, TRUE, 1L
  ),
  !.rc_step2_stream_cross_target_reuse(
    "highs_persistent_cpp", 2L, FALSE, 0L
  )
)

run_streaming_worker_check <- function(n_targets = 17L) {
  root <- tempfile("regcompass-streaming-review-")
  dir.create(root, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

  targets <- paste0("T", seq_len(n_targets))
  reactions <- c("UP", targets)
  S <- Matrix::Matrix(
    matrix(
      c(1, rep(-1, n_targets)),
      nrow = 1,
      dimnames = list("M", reactions)
    ),
    sparse = TRUE
  )
  lb <- stats::setNames(rep(0, length(reactions)), reactions)
  ub <- stats::setNames(rep(10, length(reactions)), reactions)
  units <- c("mc1", "mc2")
  primary_penalty <- matrix(
    seq(0.10, 0.90, length.out = length(reactions) * length(units)),
    nrow = length(reactions),
    dimnames = list(reactions, units)
  )
  control_penalty <- primary_penalty
  control_penalty[, "mc2"] <- rev(control_penalty[, "mc2"])
  control_identical <- c(mc1 = TRUE, mc2 = FALSE)

  row_ids <- paste0("row_", targets)
  entries <- stats::setNames(lapply(targets, function(target) {
    list(
      reaction_id = target,
      target_direction = "forward",
      cell_type = "T",
      medium_scenario = "toy"
    )
  }), row_ids)
  vmax_values <- stats::setNames(lapply(targets, function(target) {
    .rc_step2_compact_vmax_value(rc_compass_vmax_directional(
      S, lb, ub, target,
      direction = "forward", solver = "highs"
    ))
  }), row_ids)

  payload <- list(
    schema_version = "regcompass_step2_compact_payload_v1",
    model = list(
      S = S,
      lb = lb,
      ub = ub,
      target_status = "ok",
      file_checksum = "streaming-review-toy",
      cell_type = "T",
      medium_scenario = "toy"
    ),
    reactions = reactions,
    units = units,
    penalty = primary_penalty,
    penalty_evidence = .rc_step2_penalty_evidence_stats(primary_penalty),
    control_penalty = control_penalty,
    control_penalty_evidence = .rc_step2_penalty_evidence_stats(control_penalty),
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

  checkpoints <- .rc_step2_reaction_batch_worker(list(
    payload_file = payload_file,
    row_ids = row_ids,
    checkpoint_dir = checkpoint_dir
  ))
  stopifnot(
    length(checkpoints) == n_targets,
    all(file.exists(checkpoints))
  )

  expected_chunk <- .rc_step2_stream_target_chunk_size(
    length(units), n_targets
  )
  expected_chunks <- ceiling(n_targets / expected_chunk)
  stopifnot(expected_chunks >= 2L)

  for (checkpoint in checkpoints) {
    result <- readRDS(checkpoint)
    target <- entries[[result$row_id]]$reaction_id
    stopifnot(
      result$engine_metrics$stream_target_chunk_size == expected_chunk,
      result$engine_metrics$stream_target_chunk_count == expected_chunks,
      result$engine_metrics$stream_chunk_target_unit_pairs <=
        expected_chunk * length(units),
      isTRUE(result$engine_metrics$shared_model_batch_engine),
      all(result$diagnostics$step2_solver_reused_across_targets),
      !is.null(result$control),
      identical(
        as.logical(result$control$reused_from_primary_by_unit),
        unname(control_identical)
      )
    )

    for (unit in units) {
      reference <- rc_compass_two_step_lp_directional(
        S, lb, ub, target, primary_penalty[, unit],
        target_direction = "forward", omega = 0.95, solver = "highs"
      )
      stopifnot(
        identical(result$feasible[[unit]], reference$feasible),
        isTRUE(all.equal(
          result$penalty[[unit]], reference$penalty, tolerance = 1e-9
        ))
      )
      if (identical(unit, "mc1")) {
        stopifnot(isTRUE(all.equal(
          result$control$penalty[[unit]],
          result$penalty[[unit]], tolerance = 0
        )))
      } else {
        control_reference <- rc_compass_two_step_lp_directional(
          S, lb, ub, target, control_penalty[, unit],
          target_direction = "forward", omega = 0.95, solver = "highs"
        )
        stopifnot(isTRUE(all.equal(
          result$control$penalty[[unit]],
          control_reference$penalty, tolerance = 1e-9
        )))
      }
    }
  }

  # A true singleton worker may use a persistent HiGHS object for repeated
  # units, but it did not reuse that solver across targets and must say so.
  singleton_root <- file.path(root, "singleton")
  dir.create(singleton_root)
  singleton_payload <- payload
  singleton_payload$entries <- payload$entries[row_ids[[1L]]]
  singleton_payload$vmax <- payload$vmax[row_ids[[1L]]]
  singleton_file <- file.path(singleton_root, "payload.rds")
  singleton_checkpoints <- file.path(singleton_root, "checkpoints")
  dir.create(singleton_checkpoints)
  saveRDS(singleton_payload, singleton_file)
  one <- .rc_step2_reaction_batch_worker(list(
    payload_file = singleton_file,
    row_ids = row_ids[[1L]],
    checkpoint_dir = singleton_checkpoints
  ))
  singleton <- readRDS(one[[1L]])
  stopifnot(
    !isTRUE(singleton$engine_metrics$shared_model_batch_engine),
    !any(singleton$diagnostics$step2_solver_reused_across_targets)
  )

  invisible(TRUE)
}

run_streaming_worker_check()
cat("Layer 2 bounded-memory streaming/Codex review checks passed.\n")
