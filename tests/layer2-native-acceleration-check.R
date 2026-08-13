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
  metrics <- .rc_compass_step2_engine_metrics(engine)
  stopifnot(
    identical(metrics$engine, "highs_persistent_cpp"),
    metrics$n_solves == length(penalty_sets),
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
      ))
    )
  }
  invisible(result)
}

run_persistent_step2_check()
run_compact_step2_worker_check()
cat("Layer 2 native acceleration regression passed.\n")
