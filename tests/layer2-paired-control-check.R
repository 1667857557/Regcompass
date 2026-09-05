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
lb <- c(UP = -10, TARGET = -10)
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

row_key <- function(direction) {
  paste0("reaction=TARGET::direction=", direction, "::medium=base")
}

mk_payload <- function(file, primary, control = NULL, full = TRUE,
                       direction = "forward", model_status = "ok") {
  row_id <- row_key(direction)
  vmax <- list()
  vmax[[row_id]] <- list(
    feasible = TRUE, vmax = 10, status = "optimal", flux = numeric()
  )
  entry <- if (full) {
    list(
      reaction_id = "TARGET", target_direction = direction,
      medium_scenario = "base", condition = "all"
    )
  } else {
    list(
      reaction_id = "TARGET", target_direction = direction,
      medium_scenario = "base", cell_type = "T"
    )
  }
  exact_mask <- if (is.null(control)) {
    FALSE
  } else {
    value <- vapply(
      seq_len(ncol(primary)),
      function(i) identical(
        control[, i, drop = TRUE], primary[, i, drop = TRUE]
      ),
      logical(1)
    )
    names(value) <- colnames(primary)
    value
  }
  saveRDS(list(
    schema_version = if (full) {
      "regcompass_full_gem_step2_compact_payload_v1"
    } else {
      "regcompass_step2_compact_payload_v1"
    },
    model = list(
      S = S, lb = lb, ub = ub, target_status = model_status,
      file_checksum = "toy", cell_type = "T",
      medium_scenario = "base", condition = "all"
    ),
    reactions = colnames(S), units = units,
    penalty = primary,
    penalty_evidence = .rc_step2_penalty_evidence_stats(primary),
    control_penalty = control,
    control_penalty_evidence = if (is.null(control)) NULL else
      .rc_step2_penalty_evidence_stats(control),
    control_identical = exact_mask,
    entries = setNames(list(entry), row_id),
    vmax = vmax, omega = 0.95, solver = "highs", flux_threshold = 1e-8
  ), file)
}

run_worker <- function(worker, primary, control = NULL, full = TRUE,
                       direction = "forward", model_status = "ok") {
  root <- tempfile("paired-control-")
  dir.create(root)
  payload <- file.path(root, "payload.rds")
  mk_payload(payload, primary, control, full, direction, model_status)
  row_id <- row_key(direction)
  files <- worker(list(
    payload_file = payload, row_ids = row_id, checkpoint_dir = root
  ))
  out <- readRDS(files[[1L]])
  unlink(root, recursive = TRUE, force = TRUE)
  out
}

score_one <- function(result, control = FALSE) {
  row_id <- result$row_id
  value <- if (control) result$control else result
  rc_compass_score_from_penalty(
    matrix(value$penalty, nrow = 1,
           dimnames = list(row_id, names(value$penalty))),
    matrix(value$feasible, nrow = 1,
           dimnames = list(row_id, names(value$feasible)))
  )
}

check_engine_state <- function(direction) {
  vmax_result <- list(feasible = TRUE, vmax = 10, status = "optimal")
  prepared <- .rc_compass_step2_prepare(
    S, lb, ub, "TARGET", vmax_result,
    target_direction = direction, omega = 0.95, flux_threshold = 1e-8
  )
  primary_evidence <- .rc_step2_penalty_evidence_stats(primary_penalty)
  control_evidence <- .rc_step2_penalty_evidence_stats(control_penalty)
  engine <- .rc_compass_step2_new_engine(
    prepared$template, "highs", persistent_required = TRUE
  )
  on.exit(.rc_compass_step2_release_engine(engine), add = TRUE)
  first <- .rc_compass_step2_route_solve(
    prepared, engine, primary_penalty, primary_evidence, 2L, units
  )
  second <- .rc_compass_step2_route_solve(
    prepared, first$engine, control_penalty, control_evidence, 2L, units
  )
  third <- .rc_compass_step2_route_solve(
    prepared, second$engine, primary_penalty, primary_evidence, 2L, units
  )
  fresh_engine <- .rc_compass_step2_new_engine(
    prepared$template, "highs", persistent_required = TRUE
  )
  fresh <- .rc_compass_step2_route_solve(
    prepared, fresh_engine, primary_penalty, primary_evidence, 2L, units
  )
  .rc_compass_step2_release_engine(fresh$engine)
  stopifnot(
    isTRUE(all.equal(third$penalty, fresh$penalty, tolerance = 1e-12)),
    identical(third$feasible, fresh$feasible),
    identical(third$evaluated, fresh$evaluated)
  )
  engine <<- third$engine
}

check_worker <- function(worker, full, direction) {
  primary_legacy <- run_worker(
    worker, primary_penalty, full = full, direction = direction
  )
  control_legacy <- run_worker(
    worker, control_penalty, full = full, direction = direction
  )
  paired <- run_worker(
    worker, primary_penalty, control_penalty,
    full = full, direction = direction
  )
  stopifnot(
    isTRUE(all.equal(paired$penalty, primary_legacy$penalty,
                     tolerance = 1e-12)),
    identical(paired$feasible, primary_legacy$feasible),
    identical(paired$evaluated, primary_legacy$evaluated),
    identical(score_one(paired), score_one(primary_legacy)),
    isTRUE(all.equal(paired$control$penalty, control_legacy$penalty,
                     tolerance = 1e-12)),
    identical(paired$control$feasible, control_legacy$feasible),
    identical(paired$control$evaluated, control_legacy$evaluated),
    identical(score_one(paired, TRUE), score_one(control_legacy)),
    paired$engine_metrics$n_solves == length(units),
    paired$control$engine_metrics$n_solves == length(units),
    paired$control$engine_metrics$n_fallback == 0L,
    !isTRUE(paired$control$reused_from_primary),
    identical(
      unname(paired$control$reused_from_primary_by_unit),
      rep(FALSE, length(units))
    ),
    isTRUE(paired$control$shared_target_engine),
    identical(unname(paired$diagnostics$objective_value),
              unname(paired$penalty)),
    identical(unname(paired$control$diagnostics$objective_value),
              unname(paired$control$penalty))
  )

  reused <- run_worker(
    worker, primary_penalty, primary_penalty,
    full = full, direction = direction
  )
  stopifnot(
    identical(reused$penalty, reused$control$penalty),
    identical(reused$feasible, reused$control$feasible),
    identical(reused$evaluated, reused$control$evaluated),
    isTRUE(reused$control$reused_from_primary),
    reused$control$engine_metrics$n_solves == 0L,
    all(reused$control$reused_from_primary_by_unit),
    all(reused$control$diagnostics$objective_identical_to_primary),
    isTRUE(reused$control$shared_target_engine)
  )

  # The TARGET penalty is unchanged for u1, but a non-target objective
  # coefficient changes by one machine epsilon. This must solve u1 and may
  # reuse only u2/u3; target-level evidence is therefore insufficient.
  partial <- primary_penalty
  partial["UP", "u1"] <- partial["UP", "u1"] + .Machine$double.eps
  stopifnot(
    identical(partial["TARGET", "u1"], primary_penalty["TARGET", "u1"]),
    !identical(partial[, "u1"], primary_penalty[, "u1"])
  )
  partial_legacy <- run_worker(
    worker, partial, full = full, direction = direction
  )
  partial_paired <- run_worker(
    worker, primary_penalty, partial,
    full = full, direction = direction
  )
  expected_reuse <- c(u1 = FALSE, u2 = TRUE, u3 = TRUE)
  stopifnot(
    isTRUE(all.equal(
      partial_paired$control$penalty, partial_legacy$penalty,
      tolerance = 1e-12
    )),
    identical(partial_paired$control$feasible, partial_legacy$feasible),
    identical(partial_paired$control$evaluated, partial_legacy$evaluated),
    identical(
      partial_paired$control$reused_from_primary_by_unit,
      expected_reuse
    ),
    partial_paired$control$engine_metrics$n_solves == 1L,
    sum(partial_paired$control$diagnostics$objective_identical_to_primary) == 2L,
    partial_paired$control$diagnostics$solver_backend[[1L]] !=
      "reused_identical_primary_objective",
    all(
      partial_paired$control$diagnostics$solver_backend[2:3] ==
        "reused_identical_primary_objective"
    ),
    isTRUE(partial_paired$control$shared_target_engine)
  )
}

# Structural target provenance must not be rewritten from Step 2 solver
# feasibility. A deliberately distinctive model status must survive unchanged.
status_check <- run_worker(
  .rc_step2_reaction_batch_worker, primary_penalty,
  full = FALSE, model_status = "structural_status_preserved"
)
stopifnot(all(
  status_check$diagnostics$target_status == "structural_status_preserved"
))

for (direction in c("forward", "reverse")) {
  check_worker(.rc_full_gem_step2_reaction_batch_worker, TRUE, direction)
  check_worker(.rc_step2_reaction_batch_worker, FALSE, direction)
  check_engine_state(direction)
}

stage <- paste(readLines("R/step_layer2.R", warn = FALSE), collapse = "\n")
stopifnot(
  grepl("control_layer1 = control_layer1", stage, fixed = TRUE),
  grepl("paired_step2_dispatch = TRUE", stage, fixed = TRUE),
  grepl("independent_solver_streams = FALSE", stage, fixed = TRUE),
  grepl("independent_lp_solves_on_shared_target_engine = TRUE", stage,
        fixed = TRUE),
  !grepl("rna_only <- run_control(", stage, fixed = TRUE),
  !grepl("rna_only <- answer", stage, fixed = TRUE),
  grepl("answer$comparison_paths <- NULL", stage, fixed = TRUE)
)

# The paired progress wrapper must never emit phase 8 before phase 6. Use the
# paired-dispatch detail string so the legacy control wrapper's phase 8 event
# does not affect this ordering assertion.
phase6 <- regexpr(
  "primary structural models and directional scores assembled",
  stage, fixed = TRUE
)[[1L]]
phase8 <- regexpr(
  "RNA-only Step 2 completed in the same target dispatch with shared Vmax",
  stage, fixed = TRUE
)[[1L]]
stopifnot(phase6 > 0L, phase8 > 0L, phase6 < phase8)

# Reuse summaries are row-keyed in both engines, so model-scoped checkpoint
# ordering cannot scramble target-level diagnostics.
full_text <- paste(readLines("R/microcompass_engine.R", warn = FALSE),
                   collapse = "\n")
cell_text <- paste(readLines(
  "R/celltype_microcompass_reaction_parallel.R", warn = FALSE
), collapse = "\n")
for (text in list(full_text, cell_text)) {
  stopifnot(
    grepl("control_reused <- setNames(logical(length(row_ids)), row_ids)",
          text, fixed = TRUE),
    grepl("control_reused[[row_id]] <- isTRUE(all(reuse_mask))",
          text, fixed = TRUE),
    grepl("rna_control_model_identical_reuse = control_reused[row_ids]",
          text, fixed = TRUE),
    grepl("one persistent HiGHS engine per", text, fixed = TRUE)
  )
}

cat("paired Layer 2 exact-equivalence and per-unit reuse checks passed\n")
