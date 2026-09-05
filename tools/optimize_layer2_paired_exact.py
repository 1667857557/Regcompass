from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return (ROOT / path).read_text()


def write(path, text):
    (ROOT / path).write_text(text)


def replace_once(text, old, new, label):
    n = text.count(old)
    if n != 1:
        raise RuntimeError(f"{label}: expected 1 match, found {n}")
    return text.replace(old, new, 1)


def replace_regex(text, pattern, repl, label):
    out, n = re.subn(pattern, repl, text, count=1, flags=re.S)
    if n != 1:
        raise RuntimeError(f"{label}: expected 1 regex match, found {n}")
    return out


# Add exact per-unit reuse to the shared route solver. Reuse is allowed only
# when the entire model-wide objective vector is bitwise-identical for a unit.
path = "R/microcompass_vmax_cache.R"
text = read(path)
if "reused_identical_primary_objective" not in text or "reuse_mask = NULL" not in text:
    route_pattern = r'''\.rc_compass_step2_route_solve <- function\(.*?\n\}\n\n(?=\.rc_compass_step2_result <- function\()'''
    route_replacement = r'''.rc_compass_step2_engine_metrics_delta <- function(after, before) {
  fields <- c("n_solves", "n_objective_updates", "n_fallback")
  value <- lapply(fields, function(field) {
    max(
      0L,
      as.integer(after[[field]] %||% 0L) -
        as.integer(before[[field]] %||% 0L)
    )
  })
  names(value) <- fields
  c(list(engine = as.character(after$engine %||% "not_run")), value)
}

.rc_compass_step2_route_solve <- function(
    prepared, engine, penalty_matrix, evidence, target_index, units,
    reuse_mask = NULL, reuse_result = NULL) {
  penalty_matrix <- as.matrix(penalty_matrix)
  units <- as.character(units)
  n_units <- length(units)
  if (!identical(colnames(penalty_matrix), units) ||
      target_index < 1L || target_index > nrow(penalty_matrix)) {
    stop("Step 2 route penalties are not aligned to the prepared target.",
         call. = FALSE)
  }
  if (length(evidence$all_finite) != n_units ||
      length(evidence$fraction) != n_units ||
      length(evidence$unavailable) != n_units) {
    stop("Step 2 route penalty evidence summary is malformed.",
         call. = FALSE)
  }
  reuse_mask <- if (is.null(reuse_mask)) {
    rep(FALSE, n_units)
  } else {
    as.logical(reuse_mask)
  }
  if (length(reuse_mask) != n_units || anyNA(reuse_mask)) {
    stop("Step 2 exact-reuse mask is not aligned to route units.",
         call. = FALSE)
  }
  if (any(reuse_mask)) {
    required <- c(
      "penalty", "vmax", "feasible", "evaluated", "solver_status",
      "solver_backend", "step1_status", "step2_status", "target_available"
    )
    if (!is.list(reuse_result) || !all(required %in% names(reuse_result)) ||
        any(vapply(reuse_result[required], length, integer(1)) != n_units)) {
      stop("Step 2 exact-reuse source is malformed.", call. = FALSE)
    }
  }

  task_penalty <- rep(NA_real_, n_units)
  task_vmax <- rep(NA_real_, n_units)
  task_feasible <- task_evaluated <- rep(FALSE, n_units)
  names(task_penalty) <- names(task_vmax) <-
    names(task_feasible) <- names(task_evaluated) <- units
  solver_status <- step1_status <- step2_status <-
    solver_backend <- rep(NA_character_, n_units)
  target_available <- is.finite(penalty_matrix[target_index, ])

  for (i in seq_along(units)) {
    if (isTRUE(reuse_mask[[i]])) {
      task_penalty[[i]] <- reuse_result$penalty[[i]]
      task_vmax[[i]] <- reuse_result$vmax[[i]]
      task_feasible[[i]] <- reuse_result$feasible[[i]]
      task_evaluated[[i]] <- reuse_result$evaluated[[i]]
      solver_status[[i]] <- reuse_result$solver_status[[i]]
      solver_backend[[i]] <- "reused_identical_primary_objective"
      step1_status[[i]] <- reuse_result$step1_status[[i]]
      step2_status[[i]] <- reuse_result$step2_status[[i]]
      next
    }

    unit_penalty <- penalty_matrix[, i]
    solver_penalty <- if (isTRUE(evidence$all_finite[[i]])) {
      unit_penalty
    } else {
      value <- unit_penalty
      value[!is.finite(value)] <- 0
      value
    }

    if (isTRUE(prepared$runnable)) {
      solved <- .rc_compass_step2_engine_solve(
        engine, solver_penalty,
        return_solution = FALSE,
        trusted_aligned = TRUE
      )
      engine <- solved$engine
      fit <- .rc_compass_step2_result(
        prepared$template, solved$answer,
        require_solution = FALSE
      )
    } else {
      fit <- prepared$result
      solved <- NULL
    }

    task_penalty[[i]] <- if (target_available[[i]]) fit$penalty else NA_real_
    task_vmax[[i]] <- fit$vmax
    task_feasible[[i]] <- isTRUE(fit$feasible)
    task_evaluated[[i]] <- isTRUE(fit$feasible) && target_available[[i]]
    solver_status[[i]] <- as.character(fit$solver_status)
    solver_backend[[i]] <- as.character(fit$solver_backend %||% "unknown")
    step1_status[[i]] <- as.character(fit$step1_status)
    step2_status[[i]] <- as.character(fit$step2_status)
    rm(unit_penalty, solver_penalty, fit, solved)
  }

  list(
    engine = engine,
    penalty = task_penalty,
    vmax = task_vmax,
    feasible = task_feasible,
    evaluated = task_evaluated,
    solver_status = solver_status,
    solver_backend = solver_backend,
    step1_status = step1_status,
    step2_status = step2_status,
    target_available = target_available
  )
}

'''
    text = replace_regex(text, route_pattern, route_replacement,
                         "replace exact-reuse route solver")
    write(path, text)


payload_old = '''    control_identical <- identical(control_penalty_matrix, penalty_matrix)'''
payload_new = '''    control_identical <- vapply(
      seq_len(ncol(penalty_matrix)),
      function(i) identical(
        control_penalty_matrix[, i, drop = TRUE],
        penalty_matrix[, i, drop = TRUE]
      ),
      logical(1)
    )
    names(control_identical) <- colnames(penalty_matrix)'''

engine_lifetime_old = '''  checkpoint_files <- character(length(row_ids))
  step2_engine <- control_step2_engine <- NULL
  on.exit({
    .rc_compass_step2_release_engine(step2_engine)
    .rc_compass_step2_release_engine(control_step2_engine)
    rm(model, payload)
    invisible(gc(verbose = FALSE, full = TRUE))
  }, add = TRUE)'''
engine_lifetime_new = '''  checkpoint_files <- character(length(row_ids))
  step2_engine <- NULL
  on.exit({
    .rc_compass_step2_release_engine(step2_engine)
    rm(model, payload)
    invisible(gc(verbose = FALSE, full = TRUE))
  }, add = TRUE)'''

control_evidence_old = '''  control_evidence <- if (paired_control) {
    payload$control_penalty_evidence %||%
      .rc_step2_penalty_evidence_stats(payload$control_penalty)
  } else {
    NULL
  }

  for (j in seq_along(row_ids)) {'''
control_evidence_new = '''  control_evidence <- if (paired_control) {
    payload$control_penalty_evidence %||%
      .rc_step2_penalty_evidence_stats(payload$control_penalty)
  } else {
    NULL
  }
  control_reuse_mask <- if (paired_control) {
    value <- as.logical(payload$control_identical)
    if (length(value) != length(units) || anyNA(value)) {
      stop("Paired RNA-control exact-reuse mask is not aligned to units.",
           call. = FALSE)
    }
    names(value) <- units
    value
  } else {
    setNames(rep(FALSE, length(units)), units)
  }

  for (j in seq_along(row_ids)) {'''

control_block_old = '''    control <- NULL
    control_metrics <- list(
      engine = "not_run", n_solves = 0L,
      n_objective_updates = 0L, n_fallback = 0L
    )
    reused_control <- FALSE
    if (paired_control) {
      if (isTRUE(payload$control_identical)) {
        control <- primary
        control$engine <- NULL
        reused_control <- TRUE
      } else {
        control_step2_engine <- new_engine()
        control <- .rc_compass_step2_route_solve(
          prepared, control_step2_engine, payload$control_penalty,
          control_evidence, target_index, units
        )
        control_step2_engine <- control$engine
        control_metrics <-
          .rc_compass_step2_engine_metrics(control_step2_engine)
      }
    }'''
control_block_new = '''    control <- NULL
    control_metrics <- list(
      engine = "not_run", n_solves = 0L,
      n_objective_updates = 0L, n_fallback = 0L
    )
    reused_control <- control_reuse_mask
    if (paired_control) {
      control <- .rc_compass_step2_route_solve(
        prepared, step2_engine, payload$control_penalty,
        control_evidence, target_index, units,
        reuse_mask = reused_control,
        reuse_result = primary
      )
      step2_engine <- control$engine
      control_metrics <- .rc_compass_step2_engine_metrics_delta(
        .rc_compass_step2_engine_metrics(step2_engine),
        primary_metrics
      )
    }'''

control_diag_old = '''        solver_backend = if (reused_control) {
          rep("reused_identical_primary_objective", length(units))
        } else {
          control$solver_backend
        },
        objective_identical_to_primary = rep(reused_control, length(units)),'''
control_diag_new = '''        solver_backend = control$solver_backend,
        objective_identical_to_primary = unname(reused_control),'''

checkpoint_control_old = '''        engine_metrics = control_metrics,
        reused_from_primary = reused_control'''
checkpoint_control_new = '''        engine_metrics = control_metrics,
        reused_from_primary = isTRUE(all(reused_control)),
        reused_from_primary_by_unit = reused_control,
        shared_target_engine = TRUE'''

release_old = '''    .rc_compass_step2_release_engine(step2_engine)
    .rc_compass_step2_release_engine(control_step2_engine)
    step2_engine <- control_step2_engine <- NULL'''
release_new = '''    .rc_compass_step2_release_engine(step2_engine)
    step2_engine <- NULL'''

for path in ["R/microcompass_engine.R", "R/celltype_microcompass_reaction_parallel.R"]:
    text = read(path)
    if "reused_from_primary_by_unit" in text:
        continue
    text = replace_once(text, payload_old, payload_new,
                        f"{path} per-unit objective identity")
    text = replace_once(text, engine_lifetime_old, engine_lifetime_new,
                        f"{path} one target engine lifetime")
    text = replace_once(text, control_evidence_old, control_evidence_new,
                        f"{path} exact reuse mask validation")
    text = replace_once(text, control_block_old, control_block_new,
                        f"{path} shared target engine control route")
    text = replace_once(text, control_diag_old, control_diag_new,
                        f"{path} per-unit control diagnostics")
    text = replace_once(text, checkpoint_control_old, checkpoint_control_new,
                        f"{path} exact reuse checkpoint metadata")
    text = replace_once(text, release_old, release_new,
                        f"{path} target engine release")

    # Report exact reuse at the unit level while preserving the existing
    # target-level all-units-reused flag for downstream compatibility.
    aggregate_old = '''      control_diagnostics[[i]] <- ctrl$diagnostics
      control_reused[[i]] <- isTRUE(ctrl$reused_from_primary)
      control_metrics <- ctrl$engine_metrics %||% list()'''
    aggregate_new = '''      control_diagnostics[[i]] <- ctrl$diagnostics
      reuse_mask <- as.logical(
        ctrl$reused_from_primary_by_unit %||%
          rep(isTRUE(ctrl$reused_from_primary), length(result$units))
      )
      if (length(reuse_mask) != length(result$units) || anyNA(reuse_mask)) {
        stop("A paired RNA-control checkpoint has a malformed reuse mask.",
             call. = FALSE)
      }
      control_reused[[i]] <- isTRUE(all(reuse_mask))
      control_metrics <- ctrl$engine_metrics %||% list()'''
    text = replace_once(text, aggregate_old, aggregate_new,
                        f"{path} aggregate unit reuse")
    metric_old = '''        n_fallback = as.integer(control_metrics$n_fallback %||% 0L),
        reused_from_primary = isTRUE(ctrl$reused_from_primary),
        stringsAsFactors = FALSE'''
    metric_new = '''        n_fallback = as.integer(control_metrics$n_fallback %||% 0L),
        reused_from_primary = isTRUE(all(reuse_mask)),
        n_reused_from_primary = as.integer(sum(reuse_mask)),
        reuse_fraction = mean(reuse_mask),
        shared_target_engine = isTRUE(ctrl$shared_target_engine),
        stringsAsFactors = FALSE'''
    text = replace_once(text, metric_old, metric_new,
                        f"{path} exact reuse metrics")
    write(path, text)


# Replace the focused regression with a stricter forward/reverse, legacy-vs-
# paired, per-unit exact-reuse and shared-engine state-leakage regression.
test = r'''options(warn = 2)
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
                       direction = "forward") {
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
    control_identical = exact_mask,
    entries = setNames(list(entry), row_id),
    vmax = vmax, omega = 0.95, solver = "highs", flux_threshold = 1e-8
  ), file)
}

run_worker <- function(worker, primary, control = NULL, full = TRUE,
                       direction = "forward") {
  root <- tempfile("paired-control-")
  dir.create(root)
  payload <- file.path(root, "payload.rds")
  mk_payload(payload, primary, control, full, direction)
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

for (direction in c("forward", "reverse")) {
  check_worker(.rc_full_gem_step2_reaction_batch_worker, TRUE, direction)
  check_worker(.rc_step2_reaction_batch_worker, FALSE, direction)
  check_engine_state(direction)
}

stage <- paste(readLines("R/step_layer2.R", warn = FALSE), collapse = "\n")
stopifnot(
  grepl("control_layer1 = control_layer1", stage, fixed = TRUE),
  grepl("paired_step2_dispatch = TRUE", stage, fixed = TRUE),
  !grepl("rna_only <- run_control(", stage, fixed = TRUE),
  !grepl("rna_only <- answer", stage, fixed = TRUE),
  grepl("answer$comparison_paths <- NULL", stage, fixed = TRUE)
)

cat("paired Layer 2 exact-equivalence and per-unit reuse checks passed\n")
'''
write("tests/layer2-paired-control-check.R", test)

print("applied exact per-unit paired-control optimization")
