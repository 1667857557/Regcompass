`%||%` <- function(x, y) if (is.null(x)) y else x

library(Matrix)
library(Rcpp)

Rcpp::sourceCpp("src/layer2_native.cpp")
.rc_corda2_scan_flux_cpp <- rc_corda2_scan_flux_cpp

source("R/00_utils.R", local = FALSE)

rc_align_bound <- function(x, reactions, default, name) {
  if (is.null(x)) {
    value <- rep(default, length(reactions))
    names(value) <- reactions
    return(value)
  }
  if (!is.null(names(x))) {
    missing <- setdiff(reactions, names(x))
    if (length(missing)) {
      stop(name, " is missing reactions: ", paste(missing, collapse = ", "))
    }
    x <- x[reactions]
  }
  x <- as.numeric(x)
  if (length(x) != length(reactions) || anyNA(x)) {
    stop(name, " does not align with reactions")
  }
  names(x) <- reactions
  x
}

source("R/lp_solver.R", local = FALSE)
source("R/microcompass_vmax_cache.R", local = FALSE)
source("R/celltype_microcompass_reaction_parallel.R", local = FALSE)
source("R/layer2_step2_model_batch.R", local = FALSE)

.rc_atomic_save_rds <- function(x, file) {
  saveRDS(x, file)
  invisible(file)
}
.rc_bind_frames_fill <- function(x) {
  if (!length(x)) return(data.frame())
  do.call(rbind, x)
}
.rc_microcompass_object_checksum <- function(x) {
  raw <- serialize(x, NULL, version = 2L)
  paste0(length(raw), "-", sum(as.integer(raw)))
}

flux <- c(0, 2, 3, 1e-9, 4, NA_real_)
class_code <- c(1L, 2L, 3L, 4L, 2L, 3L)
track_code <- c(2L, 3L)
scan <- .rc_corda2_scan_flux_cpp(flux, class_code, track_code, 1e-8)
expected_active <- which(is.finite(flux) & flux > 1e-8)
expected_used <- expected_active[class_code[expected_active] %in% track_code]
stopifnot(
  identical(as.integer(scan$active), as.integer(expected_active)),
  identical(as.integer(scan$used), as.integer(expected_used))
)

row_ids <- paste0("row_", seq_len(12L))
model_keys <- stats::setNames(
  rep("/tmp/one-shared-model.rds", length(row_ids)), row_ids
)
tasks <- .rc_microcompass_vmax_tasks(model_keys, workers = 4L)
observed <- unlist(lapply(tasks, `[[`, "row_ids"), use.names = FALSE)
stopifnot(
  length(tasks) == 4L,
  setequal(observed, row_ids),
  length(observed) == length(unique(observed))
)

rows_a <- paste0("A", seq_len(3223L))
rows_b <- paste0("B", seq_len(4484L))
step2_model_keys <- c(
  stats::setNames(rep("model_A", length(rows_a)), rows_a),
  stats::setNames(rep("model_B", length(rows_b)), rows_b)
)
step2_tasks <- .rc_step2_model_batches(step2_model_keys, workers = 80L)
step2_observed <- unlist(
  lapply(step2_tasks, `[[`, "row_ids"), use.names = FALSE
)
step2_sizes <- vapply(
  step2_tasks, function(task) length(task$row_ids), integer(1)
)
step2_models <- vapply(step2_tasks, `[[`, character(1), "model_key")
stopifnot(
  length(step2_tasks) == 80L,
  setequal(step2_observed, names(step2_model_keys)),
  length(step2_observed) == length(unique(step2_observed)),
  max(step2_sizes) - min(step2_sizes) <= 4L,
  sum(step2_models == "model_B") > sum(step2_models == "model_A")
)

run_persistent_vmax_check <- function() {
  if (!.rc_microcompass_highs_api_available()) return(invisible(NULL))
  S <- Matrix::Matrix(
    matrix(
      c(1, -1, -1), nrow = 1,
      dimnames = list("M", c("UP", "T1", "T2"))
    ),
    sparse = TRUE
  )
  lb <- c(UP = 0, T1 = 0, T2 = 0)
  ub <- c(UP = 10, T1 = 10, T2 = 10)
  batch <- .rc_compass_vmax_batch_highs(
    S, lb, ub,
    target_reaction = c("T1", "T2"),
    direction = c("forward", "forward"),
    flux_threshold = 1e-8
  )
  stopifnot(length(batch) == 2L)
  for (i in seq_along(batch)) {
    target <- c("T1", "T2")[[i]]
    reference <- rc_compass_vmax_directional(
      S, lb, ub, target, direction = "forward", solver = "highs"
    )
    stopifnot(
      identical(batch[[i]]$feasible, reference$feasible),
      identical(batch[[i]]$status, reference$status),
      isTRUE(all.equal(batch[[i]]$vmax, reference$vmax, tolerance = 1e-10)),
      length(batch[[i]]$flux) == 0L
    )
  }
  invisible(batch)
}

run_persistent_step2_check <- function() {
  S <- Matrix::Matrix(
    matrix(
      c(1, -1), nrow = 1,
      dimnames = list("M", c("UP", "TARGET"))
    ),
    sparse = TRUE
  )
  lb <- c(UP = 0, TARGET = 0)
  ub <- c(UP = 10, TARGET = 10)
  vmax <- rc_compass_vmax_directional(
    S, lb, ub, "TARGET", direction = "forward", solver = "highs"
  )
  stopifnot(isTRUE(vmax$feasible))
  prepared <- .rc_compass_step2_prepare(
    S, lb, ub, "TARGET", vmax,
    target_direction = "forward", omega = 0.95
  )
  engine <- .rc_compass_step2_new_engine(prepared$template, "highs")
  on.exit(.rc_compass_step2_release_engine(engine), add = TRUE)
  if (!identical(engine$type, "highs_persistent_cpp")) {
    stop(
      "Persistent HiGHS Step 2 API was not activated: ",
      engine$persistent_message %||% "unknown reason"
    )
  }

  penalty_sets <- list(
    c(UP = 0.25, TARGET = 0.50),
    c(UP = 0.75, TARGET = 0.10),
    c(UP = 0.05, TARGET = 1.25)
  )
  for (penalties in penalty_sets) {
    solved <- .rc_compass_step2_engine_solve(engine, penalties)
    engine <- solved$engine
    persistent <- .rc_compass_step2_result(
      prepared$template, solved$answer
    )
    reference <- rc_compass_two_step_lp_directional(
      S, lb, ub, "TARGET", penalties,
      target_direction = "forward", omega = 0.95, solver = "highs"
    )
    stopifnot(
      identical(persistent$feasible, reference$feasible),
      isTRUE(all.equal(persistent$vmax, reference$vmax, tolerance = 1e-10)),
      isTRUE(all.equal(
        persistent$penalty, reference$penalty, tolerance = 1e-9
      ))
    )
  }
  score_only_solved <- .rc_compass_step2_engine_solve(
    engine, penalty_sets[[2L]],
    return_solution = FALSE,
    trusted_aligned = TRUE
  )
  engine <- score_only_solved$engine
  score_only <- .rc_compass_step2_result(
    prepared$template, score_only_solved$answer,
    require_solution = FALSE
  )
  score_reference <- rc_compass_two_step_lp_directional(
    S, lb, ub, "TARGET", penalty_sets[[2L]],
    target_direction = "forward", omega = 0.95, solver = "highs"
  )
  stopifnot(
    length(score_only_solved$answer$solution) == 0L,
    identical(score_only$feasible, score_reference$feasible),
    isTRUE(all.equal(
      score_only$penalty, score_reference$penalty, tolerance = 1e-9
    )),
    length(score_only$flux) == 0L
  )

  metrics <- .rc_compass_step2_engine_metrics(engine)
  stopifnot(
    identical(metrics$engine, "highs_persistent_cpp"),
    metrics$n_solves == length(penalty_sets) + 1L,
    metrics$n_objective_updates > 0L,
    metrics$n_fallback == 0L
  )
  invisible(metrics)
}

run_compact_step2_worker_check <- function() {
  root <- tempfile("regcompass-step2-worker-")
  dir.create(root, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

  S <- Matrix::Matrix(
    matrix(
      c(1, -1), nrow = 1,
      dimnames = list("M", c("UP", "TARGET"))
    ),
    sparse = TRUE
  )
  lb <- c(UP = 0, TARGET = 0)
  ub <- c(UP = 10, TARGET = 10)
  row_id <- "cell_type=T::reaction=TARGET::direction=forward::medium=toy"
  units <- c("mc1", "mc2")
  penalty_matrix <- matrix(
    c(0.25, 0.50, 0.75, 0.10),
    nrow = 2,
    dimnames = list(c("UP", "TARGET"), units)
  )
  vmax_result <- rc_compass_vmax_directional(
    S, lb, ub, "TARGET", direction = "forward", solver = "highs"
  )
  stopifnot(isTRUE(vmax_result$feasible))
  payload <- list(
    schema_version = "regcompass_step2_compact_payload_v1",
    model = list(
      S = S,
      lb = lb,
      ub = ub,
      target_status = "ok",
      file_checksum = "toy-checksum",
      cell_type = "T",
      medium_scenario = "toy"
    ),
    reactions = c("UP", "TARGET"),
    units = units,
    penalty = penalty_matrix,
    entries = stats::setNames(list(list(
      reaction_id = "TARGET",
      target_direction = "forward",
      cell_type = "T",
      medium_scenario = "toy"
    )), row_id),
    vmax = stats::setNames(list(
      .rc_step2_compact_vmax_value(vmax_result)
    ), row_id),
    omega = 0.95,
    solver = "highs",
    flux_threshold = 1e-8
  )
  payload_file <- file.path(root, "payload.rds")
  checkpoint_dir <- file.path(root, "checkpoints")
  dir.create(checkpoint_dir)
  saveRDS(payload, payload_file)

  checkpoint <- .rc_step2_reaction_batch_worker(list(
    payload_file = payload_file,
    row_ids = row_id,
    checkpoint_dir = checkpoint_dir
  ))
  stopifnot(length(checkpoint) == 1L, file.exists(checkpoint[[1L]]))
  result <- readRDS(checkpoint[[1L]])
  for (one_unit in units) {
    reference <- rc_compass_two_step_lp_directional(
      S, lb, ub, "TARGET", penalty_matrix[, one_unit],
      target_direction = "forward", omega = 0.95, solver = "highs"
    )
    stopifnot(
      identical(result$feasible[[one_unit]], reference$feasible),
      isTRUE(all.equal(
        result$vmax[[one_unit]], reference$vmax, tolerance = 1e-10
      )),
      isTRUE(all.equal(
        result$penalty[[one_unit]], reference$penalty, tolerance = 1e-9
      )),
      isTRUE(result$engine_metrics$shared_model_batch_engine),
      result$engine_metrics$batch_objective_change_events <= length(units)
    )
  }
  invisible(result)
}

run_model_batch_step2_check <- function() {
  root <- tempfile("regcompass-step2-model-batch-")
  dir.create(root, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

  S <- Matrix::Matrix(
    matrix(
      c(1, -1, -1), nrow = 1,
      dimnames = list("M", c("UP", "T1", "T2"))
    ),
    sparse = TRUE
  )
  lb <- c(UP = 0, T1 = 0, T2 = 0)
  ub <- c(UP = 10, T1 = 10, T2 = 10)
  units <- c("mc1", "mc2", "mc3")
  primary_penalty <- matrix(
    c(
      0.20, 0.50, 0.90,
      0.30, 0.70, 0.40,
      0.80, 0.15, 0.60
    ),
    nrow = 3,
    dimnames = list(c("UP", "T1", "T2"), units)
  )
  control_penalty <- primary_penalty
  control_penalty[, "mc2"] <- c(0.45, 0.25, 0.95)
  control_penalty[, "mc3"] <- c(0.10, 0.85, 0.35)
  control_identical <- c(mc1 = TRUE, mc2 = FALSE, mc3 = FALSE)
  row_ids <- c("row_T1", "row_T2")
  entries <- stats::setNames(lapply(c("T1", "T2"), function(target) {
    list(
      reaction_id = target,
      target_direction = "forward",
      cell_type = "T",
      medium_scenario = "toy"
    )
  }), row_ids)
  vmax_values <- stats::setNames(lapply(c("T1", "T2"), function(target) {
    .rc_step2_compact_vmax_value(rc_compass_vmax_directional(
      S, lb, ub, target, direction = "forward", solver = "highs"
    ))
  }), row_ids)
  payload <- list(
    schema_version = "regcompass_step2_compact_payload_v1",
    model = list(
      S = S,
      lb = lb,
      ub = ub,
      target_status = "ok",
      file_checksum = "shared-toy-model",
      cell_type = "T",
      medium_scenario = "toy"
    ),
    reactions = c("UP", "T1", "T2"),
    units = units,
    penalty = primary_penalty,
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

  checkpoints <- .rc_step2_reaction_batch_worker(list(
    payload_file = payload_file,
    row_ids = row_ids,
    checkpoint_dir = checkpoint_dir
  ))
  stopifnot(length(checkpoints) == length(row_ids), all(file.exists(checkpoints)))

  for (checkpoint in checkpoints) {
    result <- readRDS(checkpoint)
    target <- entries[[result$row_id]]$reaction_id
    stopifnot(
      result$engine_metrics$n_solves == length(units),
      isTRUE(result$engine_metrics$shared_model_batch_engine),
      result$engine_metrics$batch_objective_change_events <= length(units),
      result$engine_metrics$batch_target_switches == length(row_ids) * length(units),
      result$control$engine_metrics$n_solves == sum(!control_identical),
      result$control$engine_metrics$batch_objective_change_events <=
        sum(!control_identical),
      identical(
        result$control$diagnostics$objective_identical_to_primary,
        unname(control_identical)
      )
    )
    for (one_unit in units) {
      primary_reference <- rc_compass_two_step_lp_directional(
        S, lb, ub, target, primary_penalty[, one_unit],
        target_direction = "forward", omega = 0.95, solver = "highs"
      )
      control_reference <- rc_compass_two_step_lp_directional(
        S, lb, ub, target, control_penalty[, one_unit],
        target_direction = "forward", omega = 0.95, solver = "highs"
      )
      stopifnot(
        identical(result$feasible[[one_unit]], primary_reference$feasible),
        isTRUE(all.equal(
          result$penalty[[one_unit]], primary_reference$penalty,
          tolerance = 1e-9
        )),
        identical(
          result$control$feasible[[one_unit]], control_reference$feasible
        ),
        isTRUE(all.equal(
          result$control$penalty[[one_unit]], control_reference$penalty,
          tolerance = 1e-9
        ))
      )
    }
  }
  invisible(checkpoints)
}

run_persistent_vmax_check()
run_persistent_step2_check()
run_compact_step2_worker_check()
run_model_batch_step2_check()
cat("Layer 2 native/model-batch acceleration regression passed.\n")
