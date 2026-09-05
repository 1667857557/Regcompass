`%||%` <- function(x, y) if (is.null(x)) y else x

suppressPackageStartupMessages({
  library(Matrix)
  library(Rcpp)
  library(highs)
})

Rcpp::sourceCpp("src/layer2_native.cpp")
.rc_corda2_scan_flux_cpp <- rc_corda2_scan_flux_cpp

source("R/00_utils.R", local = FALSE)

# Standalone CI helper matching the package bound-alignment contract.
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
# This final compiler definition fixes the zero-length directional-ID edge case
# while preserving the exact PR #338 LP semantics.
source("R/layer2_step2_directional_prepare.R", local = FALSE)
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

# Preserve the existing CORDA2 native flux-scan numerical contract.  The new
# Step 2 target-batch solver never changes CORDA2 solver ownership or support
# semantics.
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

run_directional_id_edge_cases <- function() {
  S <- Matrix::Matrix(
    matrix(c(1, -1), nrow = 1,
           dimnames = list("M", c("UP", "TARGET"))),
    sparse = TRUE
  )

  # Forward-only GEM: no reverse variable may be fabricated when reverse_index
  # is integer(0).
  lb_f <- c(UP = 0, TARGET = 0)
  ub_f <- c(UP = 10, TARGET = 10)
  vmax_f <- rc_compass_vmax_directional(
    S, lb_f, ub_f, "TARGET", direction = "forward", solver = "highs"
  )
  prep_f <- .rc_compass_step2_prepare(
    S, lb_f, ub_f, "TARGET", vmax_f,
    target_direction = "forward", omega = 0.95
  )
  stopifnot(
    prep_f$runnable,
    ncol(prep_f$template$A) == length(prep_f$template$variable_id),
    ncol(prep_f$template$A) == prep_f$template$n_variables,
    all(grepl("::forward$", prep_f$template$variable_id)),
    !any(grepl("::reverse$", prep_f$template$variable_id))
  )

  # Reverse-only GEM: symmetric edge case.
  lb_r <- c(UP = -10, TARGET = -10)
  ub_r <- c(UP = 0, TARGET = 0)
  vmax_r <- rc_compass_vmax_directional(
    S, lb_r, ub_r, "TARGET", direction = "reverse", solver = "highs"
  )
  prep_r <- .rc_compass_step2_prepare(
    S, lb_r, ub_r, "TARGET", vmax_r,
    target_direction = "reverse", omega = 0.95
  )
  stopifnot(
    prep_r$runnable,
    ncol(prep_r$template$A) == length(prep_r$template$variable_id),
    ncol(prep_r$template$A) == prep_r$template$n_variables,
    all(grepl("::reverse$", prep_r$template$variable_id)),
    !any(grepl("::forward$", prep_r$template$variable_id))
  )

  invisible(TRUE)
}

run_persistent_step2_oracle_check <- function() {
  S <- Matrix::Matrix(
    matrix(c(1, -1), nrow = 1,
           dimnames = list("M", c("UP", "TARGET"))),
    sparse = TRUE
  )
  lb <- c(UP = 0, TARGET = 0)
  ub <- c(UP = 10, TARGET = 10)
  vmax <- rc_compass_vmax_directional(
    S, lb, ub, "TARGET", direction = "forward", solver = "highs"
  )
  prepared <- .rc_compass_step2_prepare(
    S, lb, ub, "TARGET", vmax,
    target_direction = "forward", omega = 0.95
  )
  engine <- .rc_compass_step2_new_engine(
    prepared$template, "highs", persistent_required = TRUE
  )
  on.exit(.rc_compass_step2_release_engine(engine), add = TRUE)
  stopifnot(identical(engine$type, "highs_persistent_cpp"))

  penalty_sets <- list(
    c(UP = 0.25, TARGET = 0.50),
    c(UP = 0.75, TARGET = 0.10),
    c(UP = 0.05, TARGET = 1.25)
  )
  for (penalties in penalty_sets) {
    solved <- .rc_compass_step2_engine_solve(
      engine, penalties, return_solution = FALSE
    )
    engine <- solved$engine
    observed <- .rc_compass_step2_result(
      prepared$template, solved$answer, require_solution = FALSE
    )
    oracle <- rc_compass_two_step_lp_directional(
      S, lb, ub, "TARGET", penalties,
      target_direction = "forward", omega = 0.95, solver = "highs"
    )
    stopifnot(
      identical(observed$feasible, oracle$feasible),
      isTRUE(all.equal(observed$vmax, oracle$vmax, tolerance = 1e-10)),
      isTRUE(all.equal(observed$penalty, oracle$penalty, tolerance = 1e-9))
    )
  }
  metrics <- .rc_compass_step2_engine_metrics(engine)
  stopifnot(
    identical(metrics$engine, "highs_persistent_cpp"),
    metrics$n_solves == length(penalty_sets),
    metrics$n_fallback == 0L
  )
  invisible(TRUE)
}

run_model_batch_step2_check <- function() {
  root <- tempfile("regcompass-step2-model-batch-")
  dir.create(root, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

  # Reversible shared model so the same persistent solver must switch between a
  # forward target and a reverse target and fully restore previous bounds.
  S <- Matrix::Matrix(
    matrix(c(1, -1, -1), nrow = 1,
           dimnames = list("M", c("UP", "T1", "T2"))),
    sparse = TRUE
  )
  lb <- c(UP = -10, T1 = -10, T2 = -10)
  ub <- c(UP = 10, T1 = 10, T2 = 10)
  reactions <- c("UP", "T1", "T2")
  units <- c("mc1", "mc2", "mc3")
  primary_penalty <- matrix(
    c(
      0.20, 0.50, 0.90,
      0.30, 0.70, 0.40,
      0.80, 0.15, 0.60
    ),
    nrow = 3,
    dimnames = list(reactions, units)
  )
  control_penalty <- primary_penalty
  control_penalty[, "mc2"] <- c(0.45, 0.25, 0.95)
  control_penalty[, "mc3"] <- c(0.10, 0.85, 0.35)
  control_identical <- c(mc1 = TRUE, mc2 = FALSE, mc3 = FALSE)

  row_ids <- c("row_T1_forward", "row_T2_reverse")
  entries <- list(
    row_T1_forward = list(
      reaction_id = "T1", target_direction = "forward",
      cell_type = "T", medium_scenario = "toy"
    ),
    row_T2_reverse = list(
      reaction_id = "T2", target_direction = "reverse",
      cell_type = "T", medium_scenario = "toy"
    )
  )
  vmax_values <- list(
    row_T1_forward = .rc_step2_compact_vmax_value(
      rc_compass_vmax_directional(
        S, lb, ub, "T1", direction = "forward", solver = "highs"
      )
    ),
    row_T2_reverse = .rc_step2_compact_vmax_value(
      rc_compass_vmax_directional(
        S, lb, ub, "T2", direction = "reverse", solver = "highs"
      )
    )
  )

  payload <- list(
    schema_version = "regcompass_step2_compact_payload_v1",
    model = list(
      S = S, lb = lb, ub = ub,
      target_status = "ok",
      file_checksum = "shared-reversible-toy-model",
      cell_type = "T", medium_scenario = "toy"
    ),
    reactions = reactions,
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
  stopifnot(
    length(checkpoints) == length(row_ids),
    all(file.exists(checkpoints))
  )

  for (checkpoint in checkpoints) {
    result <- readRDS(checkpoint)
    entry <- entries[[result$row_id]]
    stopifnot(
      result$engine_metrics$n_solves == length(units),
      isTRUE(result$engine_metrics$shared_model_batch_engine),
      result$engine_metrics$batch_objective_change_events <= length(units),
      result$engine_metrics$batch_target_switches ==
        length(row_ids) * length(units),
      result$control$engine_metrics$n_solves == sum(!control_identical),
      result$control$engine_metrics$batch_objective_change_events <=
        sum(!control_identical),
      identical(
        result$control$diagnostics$objective_identical_to_primary,
        unname(control_identical)
      )
    )

    for (one_unit in units) {
      oracle_primary <- rc_compass_two_step_lp_directional(
        S, lb, ub, entry$reaction_id, primary_penalty[, one_unit],
        target_direction = entry$target_direction,
        omega = 0.95, solver = "highs"
      )
      oracle_control <- rc_compass_two_step_lp_directional(
        S, lb, ub, entry$reaction_id, control_penalty[, one_unit],
        target_direction = entry$target_direction,
        omega = 0.95, solver = "highs"
      )
      stopifnot(
        identical(result$feasible[[one_unit]], oracle_primary$feasible),
        isTRUE(all.equal(
          result$penalty[[one_unit]], oracle_primary$penalty,
          tolerance = 1e-9
        )),
        identical(
          result$control$feasible[[one_unit]], oracle_control$feasible
        ),
        isTRUE(all.equal(
          result$control$penalty[[one_unit]], oracle_control$penalty,
          tolerance = 1e-9
        ))
      )
    }
  }

  # Existing architecture requirement: target/reaction batches remain the outer
  # parallel surface; this optimization must not serialize an 80-worker run.
  rows <- paste0("row_", seq_len(200L))
  model_keys <- stats::setNames(rep("shared_model", length(rows)), rows)
  batches <- .rc_step2_model_batches(model_keys, workers = 80L)
  assigned <- unlist(lapply(batches, `[[`, "row_ids"), use.names = FALSE)
  stopifnot(
    length(batches) == 80L,
    length(assigned) == length(rows),
    !anyDuplicated(assigned),
    setequal(assigned, rows)
  )

  invisible(TRUE)
}

run_directional_id_edge_cases()
run_persistent_step2_oracle_check()
run_model_batch_step2_check()
cat("Layer 2 exact model-batch acceleration regression passed.\n")
