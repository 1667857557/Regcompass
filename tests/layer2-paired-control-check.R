options(warn = 2)
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

source("R/layer2_parallel_runtime.R", local = FALSE)
source("R/lp_solver.R", local = FALSE)
source("R/layer2_parallel_runtime.R", local = FALSE)
source("R/microcompass_vmax_cache.R", local = FALSE)
source("R/microcompass_engine.R", local = FALSE)
source("R/celltype_microcompass_reaction_parallel.R", local = FALSE)
source("R/layer2_penalty_lp.R", local = FALSE)

if (!requireNamespace("highs", quietly = TRUE)) {
  stop("highs is required for paired-control regression")
}

S <- Matrix::Matrix(
  matrix(c(1, -1), nrow = 1,
         dimnames = list("M" = c(), c("UP", "TARGET"))),
  sparse = TRUE
)
rownames(S) <- "M"
colnames(S) <- c("UP", "TARGET")
lb <- c(UP = 0, TARGET = 0)
ub <- c(UP = 10, TARGET = 10)
units <- c("u1", "u2", "u3")
primary_penalty <- matrix(
  c(0.2, 0.7, 0.4, 0.9, 0.8, 0.3),
  nrow = 2, byrow = TRUE,
  dimnames = list(c("UP", "TARGET"), units)
)
control_penalty <- matrix(
  c(0.4, 0.6, 0.5, 0.7, 0.9, 0.2),
  nrow = 2, byrow = TRUE,
  dimnames = list(c("UP", "TARGET"), units)
)
row_id <- "reaction=TARGET::direction=forward::medium=base"
vmax <- list()
vmax[[row_id]] <- list(feasible = TRUE, vmax = 10, status = "optimal", flux = numeric())
entry_full <- list(
  reaction_id = "TARGET", target_direction = "forward",
  medium_scenario = "base", condition = "all"
)
entry_cell <- list(
  reaction_id = "TARGET", target_direction = "forward",
  medium_scenario = "base", cell_type = "T"
)

mk_payload <- function(file, primary, control = NULL, full = TRUE,
                       identical_control = FALSE) {
  saveRDS(list(
    schema_version = if (full) {
      "regcompass_full_gem_step2_compact_payload_v1"
    } else {
      "regcompass_step2_compact_payload_v1"
    },
    model = list(
      S = S, lb = lb, ub = ub, target_status = "ok",
      file_checksum = "toy", cell_type = "T",
      medium_scenario = "base", condition = "all"
    ),
    reactions = colnames(S), units = units,
    penalty = primary,
    penalty_evidence = .rc_step2_penalty_evidence_stats(primary),
    control_penalty = control,
    control_penalty_evidence = if (is.null(control)) NULL else
      .rc_step2_penalty_evidence_stats(control),
    control_identical = identical_control,
    entries = setNames(list(if (full) entry_full else entry_cell), row_id),
    vmax = vmax, omega = 0.95, solver = "highs", flux_threshold = 1e-8
  ), file)
}

run_worker <- function(worker, primary, control = NULL, full = TRUE,
                       identical_control = FALSE) {
  root <- tempfile("paired-control-")
  dir.create(root)
  payload <- file.path(root, "payload.rds")
  mk_payload(payload, primary, control, full, identical_control)
  files <- worker(list(
    payload_file = payload, row_ids = row_id, checkpoint_dir = root
  ))
  out <- readRDS(files[[1L]])
  unlink(root, recursive = TRUE, force = TRUE)
  out
}

check_worker <- function(worker, full) {
  primary_legacy <- run_worker(worker, primary_penalty, full = full)
  control_legacy <- run_worker(worker, control_penalty, full = full)
  paired <- run_worker(
    worker, primary_penalty, control_penalty, full = full
  )
  primary_score <- rc_compass_score_from_penalty(
    matrix(primary_legacy$penalty, nrow = 1,
           dimnames = list(row_id, names(primary_legacy$penalty))),
    matrix(primary_legacy$feasible, nrow = 1,
           dimnames = list(row_id, names(primary_legacy$feasible)))
  )
  paired_primary_score <- rc_compass_score_from_penalty(
    matrix(paired$penalty, nrow = 1,
           dimnames = list(row_id, names(paired$penalty))),
    matrix(paired$feasible, nrow = 1,
           dimnames = list(row_id, names(paired$feasible)))
  )
  control_score <- rc_compass_score_from_penalty(
    matrix(control_legacy$penalty, nrow = 1,
           dimnames = list(row_id, names(control_legacy$penalty))),
    matrix(control_legacy$feasible, nrow = 1,
           dimnames = list(row_id, names(control_legacy$feasible)))
  )
  paired_control_score <- rc_compass_score_from_penalty(
    matrix(paired$control$penalty, nrow = 1,
           dimnames = list(row_id, names(paired$control$penalty))),
    matrix(paired$control$feasible, nrow = 1,
           dimnames = list(row_id, names(paired$control$feasible)))
  )
  stopifnot(
    isTRUE(all.equal(paired$penalty, primary_legacy$penalty,
                     tolerance = 1e-12)),
    isTRUE(all.equal(paired$feasible, primary_legacy$feasible)),
    isTRUE(all.equal(paired$evaluated, primary_legacy$evaluated)),
    identical(paired_primary_score, primary_score),
    isTRUE(all.equal(paired$control$penalty, control_legacy$penalty,
                     tolerance = 1e-12)),
    isTRUE(all.equal(paired$control$feasible, control_legacy$feasible)),
    isTRUE(all.equal(paired$control$evaluated, control_legacy$evaluated)),
    identical(paired_control_score, control_score),
    paired$engine_metrics$n_solves == length(units),
    paired$control$engine_metrics$n_solves == length(units),
    !isTRUE(paired$control$reused_from_primary),
    identical(
      unname(paired$diagnostics$objective_value),
      unname(paired$penalty)
    ),
    identical(
      unname(paired$control$diagnostics$objective_value),
      unname(paired$control$penalty)
    )
  )

  reused <- run_worker(
    worker, primary_penalty, primary_penalty,
    full = full, identical_control = TRUE
  )
  stopifnot(
    identical(reused$penalty, reused$control$penalty),
    identical(reused$feasible, reused$control$feasible),
    identical(reused$evaluated, reused$control$evaluated),
    isTRUE(reused$control$reused_from_primary),
    reused$control$engine_metrics$n_solves == 0L
  )

  # Guard against tolerance-based reuse. A one-ULP-sized objective change must
  # still execute an independent RNA-control solver stream because the matrices
  # are not exactly identical().
  almost_control <- primary_penalty
  almost_control[[1L, 1L]] <-
    almost_control[[1L, 1L]] + .Machine$double.eps
  stopifnot(!identical(almost_control, primary_penalty))
  near <- run_worker(
    worker, primary_penalty, almost_control,
    full = full, identical_control = FALSE
  )
  stopifnot(
    !isTRUE(near$control$reused_from_primary),
    near$control$engine_metrics$n_solves == length(units)
  )
}

check_worker(.rc_full_gem_step2_reaction_batch_worker, TRUE)
check_worker(.rc_step2_reaction_batch_worker, FALSE)

stage <- paste(readLines("R/step_layer2.R", warn = FALSE), collapse = "\n")
stopifnot(
  grepl("control_layer1 = control_layer1", stage, fixed = TRUE),
  grepl("paired_step2_dispatch = TRUE", stage, fixed = TRUE),
  !grepl("rna_only <- run_control(", stage, fixed = TRUE),
  !grepl("rna_only <- answer", stage, fixed = TRUE),
  grepl("answer$comparison_paths <- NULL", stage, fixed = TRUE)
)

cat("paired Layer 2 primary/RNA-control exact-equivalence checks passed\n")