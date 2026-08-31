suppressPackageStartupMessages({
  library(Matrix)
  library(highs)
})

`%||%` <- function(x, y) if (is.null(x)) y else x
.rc_as_dgCMatrix <- function(x) methods::as(x, "dgCMatrix")
.rc_bind_frames_fill <- function(frames) {
  frames <- frames[vapply(frames, is.data.frame, logical(1))]
  if (!length(frames)) return(data.frame())
  columns <- unique(unlist(lapply(frames, colnames), use.names = FALSE))
  frames <- lapply(frames, function(x) {
    missing <- setdiff(columns, colnames(x))
    for (name in missing) x[[name]] <- NA
    x[, columns, drop = FALSE]
  })
  do.call(rbind, frames)
}

.rc_lp_status <- function(message = "", code = NA_integer_) {
  text <- tolower(paste(message, collapse = " "))
  if (grepl("infeasible", text)) return("infeasible")
  if (grepl("unbounded", text)) return("unbounded")
  if (grepl("time|limit", text)) return("time_limit")
  if (grepl("optimal", text)) return("optimal")
  if (is.finite(code) && as.integer(code) == 0L) return("optimal")
  "error"
}

rc_solve_lp <- function(obj, A, lhs, rhs, lb, ub,
                        solver = "highs", time_limit = Inf) {
  stopifnot(identical(solver, "highs"))
  answer <- highs::highs_solve(
    L = as.numeric(obj),
    lower = as.numeric(lb),
    upper = as.numeric(ub),
    A = A,
    lhs = as.numeric(lhs),
    rhs = as.numeric(rhs),
    maximum = FALSE,
    control = highs::highs_control(
      log_to_console = FALSE,
      output_flag = FALSE,
      threads = 1L,
      solver = "simplex",
      primal_feasibility_tolerance = 1e-7,
      time_limit = as.numeric(time_limit)
    )
  )
  list(
    status = .rc_lp_status(answer$status_message, answer$status),
    solution = as.numeric(answer$primal_solution),
    objective = as.numeric(answer$objective_value),
    solver_message = answer$status_message
  )
}

rc_align_bound <- function(x, rxns, default, name, allow_partial = FALSE) {
  if (is.null(x)) {
    return(stats::setNames(rep(default, length(rxns)), rxns))
  }
  if (!is.null(names(x))) x <- x[rxns]
  x <- as.numeric(x)
  if (length(x) != length(rxns) || anyNA(x)) {
    stop("invalid ", name, " bounds", call. = FALSE)
  }
  stats::setNames(x, rxns)
}

source("R/microcompass.R")
source("R/layer2_parallel_runtime.R")
source("R/microcompass_vmax_cache.R")
source("R/celltype_microcompass_reaction_parallel.R")
source("R/microcompass_engine.R")

root <- tempfile("full-gem-step2-checkpoint-ci-")
dir.create(root, recursive = TRUE)
on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

S <- Matrix::Matrix(
  matrix(
    c(1, -1, -1),
    nrow = 1,
    dimnames = list("M", c("UP", "T1", "T2"))
  ),
  sparse = TRUE
)
lb <- c(UP = 0, T1 = 0, T2 = 0)
ub <- c(UP = 10, T1 = 10, T2 = 10)
model <- list(
  S = S,
  lb = lb,
  ub = ub,
  target_status = "not_prechecked",
  closure_diagnostics = data.frame()
)
model_file <- file.path(root, "full_gem.rds")
saveRDS(model, model_file)

row_ids <- c(
  "reaction=T1::direction=forward::medium=toy::condition=all",
  "reaction=T2::direction=forward::medium=toy::condition=all"
)
model_cache <- stats::setNames(lapply(c("T1", "T2"), function(target) {
  list(
    reaction_id = target,
    target_direction = "forward",
    medium_scenario = "toy",
    condition = "all",
    file = model_file,
    build_strategy = "compass_medium_constrained_full_gem"
  )
}), row_ids)

units <- c("mc1", "mc2", "mc3")
penalty_matrix <- matrix(
  c(
    0.20, 0.50, 0.90,
    0.25, 0.80, 0.40,
    0.70, 0.10, 0.60
  ),
  nrow = 3,
  dimnames = list(c("UP", "T1", "T2"), units)
)
vmax_cache <- stats::setNames(lapply(row_ids, function(row_id) {
  list(feasible = TRUE, vmax = 10, status = "optimal", flux = numeric())
}), row_ids)

payload_dir <- file.path(root, "payload")
checkpoint_dir <- file.path(root, "checkpoint")
dir.create(payload_dir)
dir.create(checkpoint_dir)

payload_file <- .rc_full_gem_step2_model_payload(
  model_key = model_file,
  row_ids = row_ids,
  model_cache = model_cache,
  penalties = list(penalty = penalty_matrix),
  vmax_cache = vmax_cache,
  omega = 0.95,
  solver = "highs",
  flux_threshold = 1e-8,
  payload_dir = payload_dir
)
payload <- readRDS(payload_file)
stopifnot(
  identical(
    payload$schema_version,
    "regcompass_full_gem_step2_compact_payload_v1"
  ),
  identical(payload$units, units),
  isTRUE(all(vapply(payload$vmax, function(x) length(x$flux) == 0L,
                    logical(1))))
)

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
  target <- model_cache[[row_id]]$reaction_id
  stopifnot(
    identical(observed$units, units),
    identical(as.integer(observed$engine_metrics$n_solves), length(units)),
    nrow(observed$diagnostics) == length(units)
  )
  for (one_unit in units) {
    reference <- .rc_compass_step2_from_vmax_directional(
      S = S,
      lb = lb,
      ub = ub,
      target_reaction = target,
      penalties = penalty_matrix[, one_unit],
      vmax_result = vmax_cache[[row_id]],
      target_direction = "forward",
      omega = 0.95,
      solver = "highs",
      flux_threshold = 1e-8
    )
    stopifnot(
      identical(observed$feasible[[one_unit]], reference$feasible),
      isTRUE(all.equal(
        observed$vmax[[one_unit]], reference$vmax,
        tolerance = 1e-10
      )),
      isTRUE(all.equal(
        observed$penalty[[one_unit]], reference$penalty,
        tolerance = 1e-9
      ))
    )
  }
}
stopifnot(setequal(observed_rows, row_ids))

# Parallel acceleration contract: when at least 80 independent directional
# targets exist under one shared model, an 80-worker cap must remain capable of
# producing 80 reaction batches. Stability fixes must not silently serialize or
# reduce this parallel task surface.
parallel_rows <- paste0("directional_row_", seq_len(200L))
parallel_model_keys <- stats::setNames(
  rep(model_file, length(parallel_rows)),
  parallel_rows
)
parallel_batches <- .rc_step2_model_batches(
  parallel_model_keys,
  workers = 80L
)
parallel_assigned <- unlist(
  lapply(parallel_batches, `[[`, "row_ids"),
  use.names = FALSE
)
stopifnot(
  length(parallel_batches) == 80L,
  length(parallel_assigned) == length(parallel_rows),
  !anyDuplicated(parallel_assigned),
  setequal(parallel_assigned, parallel_rows)
)

cat(
  "Full-GEM direct canonical Step 2 numerical and 80-worker batching regression passed.\n"
)
