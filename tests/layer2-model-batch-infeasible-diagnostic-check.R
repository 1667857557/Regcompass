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

root <- tempfile("regcompass-step2-infeasible-diagnostic-")
dir.create(root, recursive = TRUE)

reactions <- c("UP", "T1", "T2", "BAD")
S <- Matrix::Matrix(
  matrix(
    c(1, -1, -1, -1),
    nrow = 1,
    dimnames = list("M", reactions)
  ),
  sparse = TRUE
)
lb <- c(UP = 0, T1 = 0, T2 = 0, BAD = 0)
ub <- c(UP = 10, T1 = 10, T2 = 10, BAD = 0)
units <- c("mc1", "mc2")
penalty <- matrix(
  c(
    0.20, 0.40, 0.60, 0.80,
    0.25, 0.45, 0.65, 0.85
  ),
  nrow = length(reactions),
  dimnames = list(reactions, units)
)
row_ids <- c("row_T1", "row_T2", "row_BAD")
entries <- list(
  row_T1 = list(
    reaction_id = "T1", target_direction = "forward",
    cell_type = "T", medium_scenario = "toy"
  ),
  row_T2 = list(
    reaction_id = "T2", target_direction = "forward",
    cell_type = "T", medium_scenario = "toy"
  ),
  row_BAD = list(
    reaction_id = "BAD", target_direction = "forward",
    cell_type = "T", medium_scenario = "toy"
  )
)
vmax <- list(
  row_T1 = .rc_step2_compact_vmax_value(rc_compass_vmax_directional(
    S, lb, ub, "T1", direction = "forward", solver = "highs"
  )),
  row_T2 = .rc_step2_compact_vmax_value(rc_compass_vmax_directional(
    S, lb, ub, "T2", direction = "forward", solver = "highs"
  )),
  row_BAD = list(
    feasible = FALSE,
    vmax = NA_real_,
    status = "infeasible",
    flux = numeric()
  )
)
stopifnot(
  isTRUE(vmax$row_T1$feasible),
  isTRUE(vmax$row_T2$feasible),
  !isTRUE(vmax$row_BAD$feasible)
)

payload <- list(
  schema_version = "regcompass_step2_compact_payload_v1",
  model = list(
    S = S,
    lb = lb,
    ub = ub,
    target_status = "mixed",
    file_checksum = "infeasible-diagnostic-toy",
    cell_type = "T",
    medium_scenario = "toy"
  ),
  reactions = reactions,
  units = units,
  penalty = penalty,
  penalty_evidence = .rc_step2_penalty_evidence_stats(penalty),
  control_penalty = NULL,
  entries = entries,
  vmax = vmax,
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
result <- stats::setNames(lapply(checkpoints, readRDS), row_ids)

for (row_id in c("row_T1", "row_T2")) {
  one <- result[[row_id]]
  stopifnot(
    isTRUE(one$engine_metrics$shared_model_batch_engine),
    all(one$diagnostics$step2_solver_reused_across_targets),
    as.integer(one$engine_metrics$n_solves) == length(units)
  )
}

bad <- result$row_BAD
stopifnot(
  !isTRUE(bad$engine_metrics$shared_model_batch_engine),
  !any(bad$diagnostics$step2_solver_reused_across_targets),
  as.integer(bad$engine_metrics$n_solves) == 0L,
  all(!bad$feasible),
  all(!bad$evaluated),
  all(is.na(bad$penalty))
)

unlink(root, recursive = TRUE, force = TRUE)
cat("Layer 2 infeasible-row reuse diagnostics check passed.\n")
