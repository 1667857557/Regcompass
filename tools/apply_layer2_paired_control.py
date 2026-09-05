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

# Version bump.
desc = read("DESCRIPTION")
desc = replace_once(desc, "Version: 2.4.23", "Version: 2.4.24", "DESCRIPTION version")
write("DESCRIPTION", desc)

# Shared route solver: one prepared target template, independent solver stream per route.
path = "R/microcompass_vmax_cache.R"
text = read(path)
anchor = '''.rc_compass_step2_result <- function(\n    template, answer, require_solution = TRUE) {'''
helper = r'''.rc_compass_step2_route_solve <- function(
    prepared, engine, penalty_matrix, evidence, target_index, units) {
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

  task_penalty <- rep(NA_real_, n_units)
  task_vmax <- rep(NA_real_, n_units)
  task_feasible <- task_evaluated <- rep(FALSE, n_units)
  names(task_penalty) <- names(task_vmax) <-
    names(task_feasible) <- names(task_evaluated) <- units
  solver_status <- step1_status <- step2_status <-
    solver_backend <- rep(NA_character_, n_units)
  target_available <- is.finite(penalty_matrix[target_index, ])

  for (i in seq_along(units)) {
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
text = replace_once(text, anchor, helper + anchor, "insert paired route solver")
write(path, text)

# Full-GEM payload + worker.
path = "R/microcompass_engine.R"
text = read(path)
full_payload = r'''.rc_full_gem_step2_model_payload <- function(
    model_key, row_ids, model_cache, penalties, vmax_cache,
    omega, solver, flux_threshold, payload_dir,
    control_penalties = NULL) {
  row_ids <- as.character(row_ids)
  if (!length(row_ids)) {
    stop("A full-GEM Step 2 payload requires at least one target row.",
         call. = FALSE)
  }
  first_entry <- model_cache[[row_ids[[1L]]]]
  if (is.null(first_entry)) {
    stop("A full-GEM Step 2 payload references an unknown target row.",
         call. = FALSE)
  }
  model <- .rc_load_microcompass_model(first_entry, "full_gem")
  reactions <- colnames(model$S)
  if (is.null(reactions) || anyNA(reactions) || any(!nzchar(reactions)) ||
      anyDuplicated(reactions)) {
    stop("A full-GEM Step 2 payload has invalid reaction identifiers.",
         call. = FALSE)
  }
  units <- colnames(penalties$penalty)
  if (is.null(units) || anyNA(units) || any(!nzchar(units)) ||
      anyDuplicated(units)) {
    stop("A full-GEM Step 2 payload has invalid unit identifiers.",
         call. = FALSE)
  }
  missing_penalty_rows <- setdiff(reactions, rownames(penalties$penalty))
  if (length(missing_penalty_rows)) {
    stop(
      "Full-GEM Step 2 penalties are missing model reactions: ",
      paste(utils::head(missing_penalty_rows, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  entries <- lapply(row_ids, function(row_id) {
    entry <- model_cache[[row_id]]
    if (is.null(entry) ||
        !identical(as.character(entry$file), as.character(model_key))) {
      stop("A full-GEM Step 2 payload mixes different model files.",
           call. = FALSE)
    }
    list(
      reaction_id = as.character(entry$reaction_id),
      target_direction = as.character(entry$target_direction),
      medium_scenario = as.character(entry$medium_scenario),
      condition = as.character(entry$condition %||% "all")
    )
  })
  names(entries) <- row_ids
  vmax_values <- lapply(row_ids, function(row_id) {
    value <- vmax_cache[[row_id]]
    list(
      feasible = isTRUE(value$feasible),
      vmax = as.numeric(value$vmax),
      status = as.character(value$status),
      flux = numeric()
    )
  })
  names(vmax_values) <- row_ids
  model_checksum <- if (!is.null(first_entry$file) &&
                        file.exists(first_entry$file)) {
    unname(tools::md5sum(first_entry$file)[[1L]])
  } else {
    NA_character_
  }
  penalty_matrix <- penalties$penalty[reactions, units, drop = FALSE]
  penalty_evidence <- .rc_step2_penalty_evidence_stats(penalty_matrix)
  control_penalty_matrix <- NULL
  control_penalty_evidence <- NULL
  control_identical <- FALSE
  if (!is.null(control_penalties)) {
    if (!identical(colnames(control_penalties$penalty), units) ||
        !all(reactions %in% rownames(control_penalties$penalty))) {
      stop("Full-GEM RNA-control penalties are not aligned to primary penalties.",
           call. = FALSE)
    }
    control_penalty_matrix <-
      control_penalties$penalty[reactions, units, drop = FALSE]
    control_penalty_evidence <-
      .rc_step2_penalty_evidence_stats(control_penalty_matrix)
    control_identical <- identical(control_penalty_matrix, penalty_matrix)
  }
  payload <- list(
    schema_version = "regcompass_full_gem_step2_compact_payload_v1",
    model = list(
      S = .rc_as_dgCMatrix(model$S),
      lb = model$lb,
      ub = model$ub,
      target_status = model$target_status %||% "not_prechecked",
      file_checksum = model_checksum,
      medium_scenario = as.character(first_entry$medium_scenario),
      condition = as.character(first_entry$condition %||% "all")
    ),
    reactions = reactions,
    units = units,
    penalty = penalty_matrix,
    penalty_evidence = penalty_evidence,
    control_penalty = control_penalty_matrix,
    control_penalty_evidence = control_penalty_evidence,
    control_identical = control_identical,
    entries = entries,
    vmax = vmax_values,
    omega = as.numeric(omega),
    solver = as.character(solver),
    flux_threshold = as.numeric(flux_threshold)
  )
  token <- substr(.rc_microcompass_object_checksum(list(
    file = as.character(first_entry$file %||% model_key),
    checksum = model_checksum,
    units = units,
    row_ids = row_ids
  )), 1L, 24L)
  file <- file.path(payload_dir, paste0("payload__", token, ".rds"))
  .rc_atomic_save_rds(payload, file)
  rm(
    model, payload, entries, vmax_values, penalty_matrix, penalty_evidence,
    control_penalty_matrix, control_penalty_evidence
  )
  invisible(gc(verbose = FALSE, full = FALSE))
  file
}

'''
text = replace_regex(
    text,
    r'\.rc_full_gem_step2_model_payload <- function\(.*?\n}\n\n(?=\.rc_full_gem_step2_reaction_batch_worker)',
    full_payload,
    "replace full-GEM payload"
)

full_worker = r'''.rc_full_gem_step2_reaction_batch_worker <- function(task) {
  if (!is.list(task) ||
      !all(c("payload_file", "row_ids", "checkpoint_dir") %in% names(task))) {
    stop("Malformed full-GEM Step 2 reaction-batch task.", call. = FALSE)
  }
  payload <- readRDS(task$payload_file)
  if (!is.list(payload) ||
      !identical(
        payload$schema_version,
        "regcompass_full_gem_step2_compact_payload_v1"
      )) {
    stop("Malformed full-GEM Step 2 compact payload.", call. = FALSE)
  }
  row_ids <- as.character(task$row_ids)
  if (!length(row_ids) ||
      !all(row_ids %in% names(payload$entries)) ||
      !all(row_ids %in% names(payload$vmax))) {
    stop("Full-GEM Step 2 reaction-batch rows are absent from the payload.",
         call. = FALSE)
  }
  model <- payload$model
  if (!is.list(model) || is.null(model$S) || is.null(model$lb) ||
      is.null(model$ub)) {
    stop("Full-GEM Step 2 compact payload lacks the required LP model state.",
         call. = FALSE)
  }
  if (!identical(colnames(model$S), as.character(payload$reactions))) {
    stop("Full-GEM Step 2 payload reaction order differs from its model.",
         call. = FALSE)
  }
  if (!identical(colnames(payload$penalty), as.character(payload$units)) ||
      !identical(rownames(payload$penalty), as.character(payload$reactions))) {
    stop("Full-GEM Step 2 compact payload penalties are not aligned.",
         call. = FALSE)
  }
  paired_control <- !is.null(payload$control_penalty)
  if (paired_control &&
      (!identical(dimnames(payload$control_penalty), dimnames(payload$penalty)))) {
    stop("Full-GEM paired RNA-control penalties are not aligned.", call. = FALSE)
  }

  checkpoint_files <- character(length(row_ids))
  step2_engine <- control_step2_engine <- NULL
  on.exit({
    .rc_compass_step2_release_engine(step2_engine)
    .rc_compass_step2_release_engine(control_step2_engine)
    rm(model, payload)
    invisible(gc(verbose = FALSE, full = TRUE))
  }, add = TRUE)

  units <- as.character(payload$units)
  evidence <- payload$penalty_evidence %||%
    .rc_step2_penalty_evidence_stats(payload$penalty)
  control_evidence <- if (paired_control) {
    payload$control_penalty_evidence %||%
      .rc_step2_penalty_evidence_stats(payload$control_penalty)
  } else {
    NULL
  }

  for (j in seq_along(row_ids)) {
    row_id <- row_ids[[j]]
    entry <- payload$entries[[row_id]]
    target_index <- match(entry$reaction_id, payload$reactions)
    if (is.na(target_index)) {
      stop("A full-GEM target reaction is absent from its shared model.",
           call. = FALSE)
    }
    prepared <- .rc_compass_step2_prepare(
      S = model$S,
      lb = model$lb,
      ub = model$ub,
      target_reaction = entry$reaction_id,
      vmax_result = payload$vmax[[row_id]],
      target_direction = entry$target_direction,
      omega = payload$omega,
      flux_threshold = payload$flux_threshold
    )
    new_engine <- function() {
      if (!isTRUE(prepared$runnable)) return(NULL)
      .rc_compass_step2_new_engine(
        prepared$template, payload$solver,
        persistent_required = identical(payload$solver, "highs")
      )
    }
    step2_engine <- new_engine()
    primary <- .rc_compass_step2_route_solve(
      prepared, step2_engine, payload$penalty, evidence,
      target_index, units
    )
    step2_engine <- primary$engine
    primary_metrics <- .rc_compass_step2_engine_metrics(step2_engine)

    control <- NULL
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
    }

    diagnostics <- data.frame(
      row_id = rep(row_id, length(units)),
      unit_id = units,
      module_id = rep(NA_character_, length(units)),
      reaction_id = rep(entry$reaction_id, length(units)),
      target_direction = rep(entry$target_direction, length(units)),
      medium_scenario = rep(entry$medium_scenario, length(units)),
      condition = rep(entry$condition, length(units)),
      strict_feasible = unname(primary$feasible),
      solver_status = primary$solver_status,
      solver_backend = primary$solver_backend,
      step1_status = primary$step1_status,
      step2_status = primary$step2_status,
      target_status = ifelse(
        primary$feasible, "ok", "medium_directionally_infeasible"
      ),
      objective_value = unname(primary$penalty),
      vmax = unname(primary$vmax),
      vmax_reused_from_shared_cache = rep(TRUE, length(units)),
      step2_model_reused_across_metacells = rep(TRUE, length(units)),
      target_expression_available = primary$target_available,
      objective_evidence_fraction = evidence$fraction,
      unavailable_objective_terms = evidence$unavailable,
      parallel_task = rep(
        "directional_reaction_x_all_metacells", length(units)
      ),
      stringsAsFactors = FALSE
    )
    control_diagnostics <- if (paired_control) {
      data.frame(
        row_id = rep(row_id, length(units)),
        unit_id = units,
        reaction_id = rep(entry$reaction_id, length(units)),
        target_direction = rep(entry$target_direction, length(units)),
        medium_scenario = rep(entry$medium_scenario, length(units)),
        objective_value = unname(control$penalty),
        strict_feasible = unname(control$feasible),
        solver_status = control$solver_status,
        solver_backend = if (reused_control) {
          rep("reused_identical_primary_objective", length(units))
        } else {
          control$solver_backend
        },
        objective_identical_to_primary = rep(reused_control, length(units)),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame()
    }

    token <- substr(.rc_microcompass_object_checksum(list(
      row_id = row_id,
      file_checksum = model$file_checksum,
      units = units,
      omega = payload$omega,
      solver = payload$solver,
      flux_threshold = payload$flux_threshold
    )), 1L, 24L)
    checkpoint <- file.path(
      task$checkpoint_dir, paste0("step2__", token, ".rds")
    )
    .rc_atomic_save_rds(list(
      row_id = row_id,
      units = units,
      penalty = primary$penalty,
      vmax = primary$vmax,
      feasible = primary$feasible,
      evaluated = primary$evaluated,
      diagnostics = diagnostics,
      engine_metrics = primary_metrics,
      control = if (paired_control) list(
        penalty = control$penalty,
        vmax = control$vmax,
        feasible = control$feasible,
        evaluated = control$evaluated,
        diagnostics = control_diagnostics,
        engine_metrics = control_metrics,
        reused_from_primary = reused_control
      ) else NULL
    ), checkpoint)
    checkpoint_files[[j]] <- checkpoint

    .rc_compass_step2_release_engine(step2_engine)
    .rc_compass_step2_release_engine(control_step2_engine)
    step2_engine <- control_step2_engine <- NULL
    rm(
      diagnostics, control_diagnostics, prepared, primary, control,
      primary_metrics, control_metrics
    )
    invisible(gc(verbose = FALSE, full = FALSE))
  }
  checkpoint_files
}

'''
text = replace_regex(
    text,
    r'\.rc_full_gem_step2_reaction_batch_worker <- function\(task\) \{.*?\n}\n\n(?=\.rc_run_shared_full_gem_engine_core)',
    full_worker,
    "replace full-GEM worker"
)

# Add optional control layer to full-GEM core and compute both penalty matrices once.
text = replace_once(
    text,
    "    BPPARAM = NULL,\n    model_cache_override = NULL) {\n  mode <- match.arg(mode)",
    "    BPPARAM = NULL,\n    model_cache_override = NULL,\n    control_layer1 = NULL) {\n  mode <- match.arg(mode)",
    "full-GEM core control arg"
)
text = replace_once(
    text,
    '''  matrices <- rc_layer2_unit_matrices(\n    layer1,\n    if (identical(unit, "metacell")) "metacell" else "sample_celltype",\n    sample_col, celltype_col, condition_col\n  )\n  gem <- rc_annotate_reaction_roles(gem)''',
    '''  matrices <- rc_layer2_unit_matrices(\n    layer1,\n    if (identical(unit, "metacell")) "metacell" else "sample_celltype",\n    sample_col, celltype_col, condition_col\n  )\n  control_matrices <- if (!is.null(control_layer1)) {\n    value <- rc_layer2_unit_matrices(\n      control_layer1,\n      if (identical(unit, "metacell")) "metacell" else "sample_celltype",\n      sample_col, celltype_col, condition_col\n    )\n    if (!identical(dimnames(value$reaction_expression),\n                   dimnames(matrices$reaction_expression)) ||\n        !identical(value$unit_meta, matrices$unit_meta)) {\n      stop("Paired RNA-control units do not align with primary Layer 2 units.",\n           call. = FALSE)\n    }\n    value\n  } else NULL\n  gem <- rc_annotate_reaction_roles(gem)''',
    "full-GEM control matrices"
)
text = replace_once(
    text,
    '''  penalties <- rc_compute_multiome_penalty(\n    rc_align_reaction_expression(\n      matrices$reaction_expression, all_reactions, NA_real_\n    ),\n    reaction_roles = gem$reaction_roles\n  )\n  vmax_cache <- .rc_build_microcompass_vmax_cache(''',
    '''  penalties <- rc_compute_multiome_penalty(\n    rc_align_reaction_expression(\n      matrices$reaction_expression, all_reactions, NA_real_\n    ),\n    reaction_roles = gem$reaction_roles\n  )\n  control_penalties <- if (!is.null(control_matrices)) {\n    rc_compute_multiome_penalty(\n      rc_align_reaction_expression(\n        control_matrices$reaction_expression, all_reactions, NA_real_\n      ),\n      reaction_roles = gem$reaction_roles\n    )\n  } else NULL\n  vmax_cache <- .rc_build_microcompass_vmax_cache(''',
    "full-GEM control penalties"
)
text = replace_once(
    text,
    '''  feasible <- evaluated <- matrix(\n    FALSE,\n    nrow = length(row_ids),\n    ncol = length(units),\n    dimnames = list(row_ids, units)\n  )\n\n  checkpoint_root''',
    '''  feasible <- evaluated <- matrix(\n    FALSE,\n    nrow = length(row_ids),\n    ncol = length(units),\n    dimnames = list(row_ids, units)\n  )\n  control_penalty <- if (!is.null(control_penalties)) penalty else NULL\n  control_feasible <- if (!is.null(control_penalties)) feasible else NULL\n  control_evaluated <- if (!is.null(control_penalties)) evaluated else NULL\n\n  checkpoint_root''',
    "full-GEM control output matrices"
)
text = replace_once(
    text,
    '''      flux_threshold = flux_threshold,\n      payload_dir = payload_dir\n    )''',
    '''      flux_threshold = flux_threshold,\n      payload_dir = payload_dir,\n      control_penalties = control_penalties\n    )''',
    "full-GEM payload control pass"
)
text = replace_once(
    text,
    '''  penalties$penalty <- NULL\n  rm(vmax_cache, matrices, all_reactions)''',
    '''  penalties$penalty <- NULL\n  if (!is.null(control_penalties)) control_penalties$penalty <- NULL\n  rm(vmax_cache, matrices, control_matrices, all_reactions)''',
    "full-GEM cleanup control"
)
text = replace_once(
    text,
    '''  diagnostics <- vector("list", length(checkpoint_files))\n  step2_engine_metrics <- vector("list", length(checkpoint_files))\n  observed_rows <- character(length(checkpoint_files))''',
    '''  diagnostics <- vector("list", length(checkpoint_files))\n  step2_engine_metrics <- vector("list", length(checkpoint_files))\n  control_diagnostics <- vector("list", length(checkpoint_files))\n  control_step2_engine_metrics <- vector("list", length(checkpoint_files))\n  control_reused <- logical(length(checkpoint_files))\n  observed_rows <- character(length(checkpoint_files))''',
    "full-GEM control collector init"
)
text = replace_once(
    text,
    '''    step2_engine_metrics[[i]] <- data.frame(\n      row_id = row_id,\n      engine = as.character(metrics$engine %||% "unknown"),\n      n_solves = as.integer(metrics$n_solves %||% 0L),\n      n_objective_updates = as.integer(\n        metrics$n_objective_updates %||% 0L\n      ),\n      n_fallback = as.integer(metrics$n_fallback %||% 0L),\n      stringsAsFactors = FALSE\n    )\n    rm(result, metrics)''',
    '''    step2_engine_metrics[[i]] <- data.frame(\n      row_id = row_id,\n      engine = as.character(metrics$engine %||% "unknown"),\n      n_solves = as.integer(metrics$n_solves %||% 0L),\n      n_objective_updates = as.integer(\n        metrics$n_objective_updates %||% 0L\n      ),\n      n_fallback = as.integer(metrics$n_fallback %||% 0L),\n      stringsAsFactors = FALSE\n    )\n    if (!is.null(control_penalty)) {\n      ctrl <- result$control\n      if (!is.list(ctrl) ||\n          !identical(names(ctrl$penalty), as.character(result$units))) {\n        stop("A paired full-GEM RNA-control checkpoint is malformed.",\n             call. = FALSE)\n      }\n      control_penalty[row_id, result$units] <- ctrl$penalty\n      control_feasible[row_id, result$units] <- ctrl$feasible\n      control_evaluated[row_id, result$units] <- ctrl$evaluated\n      control_diagnostics[[i]] <- ctrl$diagnostics\n      control_reused[[i]] <- isTRUE(ctrl$reused_from_primary)\n      control_metrics <- ctrl$engine_metrics %||% list()\n      control_step2_engine_metrics[[i]] <- data.frame(\n        row_id = row_id,\n        engine = as.character(control_metrics$engine %||% "not_run"),\n        n_solves = as.integer(control_metrics$n_solves %||% 0L),\n        n_objective_updates = as.integer(\n          control_metrics$n_objective_updates %||% 0L\n        ),\n        n_fallback = as.integer(control_metrics$n_fallback %||% 0L),\n        reused_from_primary = isTRUE(ctrl$reused_from_primary),\n        stringsAsFactors = FALSE\n      )\n    }\n    rm(result, metrics)''',
    "full-GEM collect control"
)
text = replace_once(
    text,
    '''  score <- rc_compass_score_from_penalty(penalty, feasible)\n  lp_diagnostics <- .rc_bind_frames_fill(diagnostics)''',
    '''  score <- rc_compass_score_from_penalty(penalty, feasible)\n  score_rna_only <- if (!is.null(control_penalty)) {\n    rc_compass_score_from_penalty(control_penalty, control_feasible)\n  } else NULL\n  lp_diagnostics <- .rc_bind_frames_fill(diagnostics)''',
    "full-GEM control score"
)
text = replace_once(
    text,
    '''    evaluated = evaluated,\n    target_direction = directions,''',
    '''    evaluated = evaluated,\n    penalty_rna_only = control_penalty,\n    score_rna_only = score_rna_only,\n    feasible_rna_only = control_feasible,\n    evaluated_rna_only = control_evaluated,\n    lp_diagnostics_rna_only = .rc_bind_frames_fill(control_diagnostics),\n    step2_engine_metrics_rna_only =\n      .rc_bind_frames_fill(control_step2_engine_metrics),\n    rna_control_model_identical_reuse = control_reused,\n    target_direction = directions,''',
    "full-GEM return control"
)
text = replace_once(
    text,
    '''      step2_solver_reuse = paste(\n        "one persistent HiGHS model per directional reaction reused across",\n        "all metacells by objective-coefficient updates"\n      ),''',
    '''      step2_solver_reuse = paste(\n        "one prepared target template with independent persistent HiGHS",\n        "streams for primary and RNA-only objectives across all metacells"\n      ),\n      paired_primary_rna_control = !is.null(control_penalty),\n      rna_control_vmax_solve_count = if (!is.null(control_penalty)) 0L else NA_integer_,\n      rna_control_identical_model_reuse = sum(control_reused),''',
    "full-GEM params control"
)

# Thread control_layer1 through internal/public microCOMPASS APIs.
text = replace_once(
    text,
    '''    BPPARAM = NULL,\n    model_cache_override = NULL) {\n  mode <- match.arg(mode)''',
    '''    BPPARAM = NULL,\n    model_cache_override = NULL,\n    control_layer1 = NULL) {\n  mode <- match.arg(mode)''',
    "microcompass engine control arg"
)
text = replace_once(
    text,
    '''      BPPARAM = BPPARAM,\n      model_cache_override = model_cache_override\n    ))''',
    '''      BPPARAM = BPPARAM,\n      model_cache_override = model_cache_override,\n      control_layer1 = control_layer1\n    ))''',
    "celltype engine control pass"
)
text = replace_once(
    text,
    '''    BPPARAM = BPPARAM,\n    model_cache_override = model_cache_override\n  )''',
    '''    BPPARAM = BPPARAM,\n    model_cache_override = model_cache_override,\n    control_layer1 = control_layer1\n  )''',
    "full engine control pass"
)
text = replace_once(
    text,
    '''    flux_threshold = 1e-8,\n    BPPARAM = NULL) {''',
    '''    flux_threshold = 1e-8,\n    BPPARAM = NULL,\n    control_layer1 = NULL) {''',
    "public microcompass control arg"
)
text = replace_once(
    text,
    '''    flux_threshold = flux_threshold,\n    BPPARAM = BPPARAM\n  )''',
    '''    flux_threshold = flux_threshold,\n    BPPARAM = BPPARAM,\n    control_layer1 = control_layer1\n  )''',
    "public microcompass control pass"
)
write(path, text)

# Cell-type payload + worker + core paired control.
path = "R/celltype_microcompass_reaction_parallel.R"
text = read(path)
cell_payload = r'''.rc_step2_model_payload <- function(
    model_key, row_ids, model_cache, unit_celltype, penalties, vmax_cache,
    omega, solver, flux_threshold, payload_dir,
    control_penalties = NULL) {
  row_ids <- as.character(row_ids)
  first_entry <- model_cache[[row_ids[[1L]]]]
  eligible <- names(unit_celltype)[unit_celltype == first_entry$cell_type]
  if (!length(eligible)) {
    stop("No Layer 1 units match cell type `", first_entry$cell_type, "`.",
         call. = FALSE)
  }
  model <- .rc_read_celltype_union_gem(
    first_entry$file, first_entry$cell_type,
    first_entry$medium_scenario, first_entry$file_checksum
  )
  reactions <- colnames(model$S)
  if (is.null(reactions) || anyNA(reactions) || any(!nzchar(reactions))) {
    stop("A cell-type union GEM has invalid reaction identifiers.",
         call. = FALSE)
  }
  entries <- lapply(row_ids, function(row_id) {
    entry <- model_cache[[row_id]]
    if (!identical(as.character(entry$file), as.character(model_key))) {
      stop("A Step 2 payload mixes different union GEM files.", call. = FALSE)
    }
    list(
      reaction_id = as.character(entry$reaction_id),
      target_direction = as.character(entry$target_direction),
      cell_type = as.character(entry$cell_type),
      medium_scenario = as.character(entry$medium_scenario)
    )
  })
  names(entries) <- row_ids
  vmax_values <- lapply(row_ids, function(row_id) {
    .rc_step2_compact_vmax_value(vmax_cache[[row_id]])
  })
  names(vmax_values) <- row_ids
  penalty_matrix <- penalties$penalty[reactions, eligible, drop = FALSE]
  penalty_evidence <- .rc_step2_penalty_evidence_stats(penalty_matrix)
  control_penalty_matrix <- NULL
  control_penalty_evidence <- NULL
  control_identical <- FALSE
  if (!is.null(control_penalties)) {
    if (!identical(colnames(control_penalties$penalty), colnames(penalties$penalty)) ||
        !all(reactions %in% rownames(control_penalties$penalty))) {
      stop("RNA-control penalties are not aligned to primary penalties.",
           call. = FALSE)
    }
    control_penalty_matrix <-
      control_penalties$penalty[reactions, eligible, drop = FALSE]
    control_penalty_evidence <-
      .rc_step2_penalty_evidence_stats(control_penalty_matrix)
    control_identical <- identical(control_penalty_matrix, penalty_matrix)
  }
  payload <- list(
    schema_version = "regcompass_step2_compact_payload_v1",
    model = list(
      S = model$S,
      lb = model$lb,
      ub = model$ub,
      target_status = model$target_status %||% NA_character_,
      file_checksum = as.character(first_entry$file_checksum),
      cell_type = as.character(first_entry$cell_type),
      medium_scenario = as.character(first_entry$medium_scenario)
    ),
    reactions = reactions,
    units = eligible,
    penalty = penalty_matrix,
    penalty_evidence = penalty_evidence,
    control_penalty = control_penalty_matrix,
    control_penalty_evidence = control_penalty_evidence,
    control_identical = control_identical,
    entries = entries,
    vmax = vmax_values,
    omega = as.numeric(omega),
    solver = as.character(solver),
    flux_threshold = as.numeric(flux_threshold)
  )
  token <- substr(.rc_microcompass_object_checksum(list(
    file = first_entry$file,
    checksum = first_entry$file_checksum,
    units = eligible,
    row_ids = row_ids
  )), 1L, 24L)
  file <- file.path(payload_dir, paste0("payload__", token, ".rds"))
  .rc_atomic_save_rds(payload, file)
  rm(
    model, payload, entries, vmax_values, penalty_matrix, penalty_evidence,
    control_penalty_matrix, control_penalty_evidence
  )
  invisible(gc(verbose = FALSE, full = FALSE))
  file
}

'''
text = replace_regex(
    text,
    r'\.rc_step2_model_payload <- function\(.*?\n}\n\n(?=\.rc_step2_reaction_batch_worker)',
    cell_payload,
    "replace celltype payload"
)
cell_worker = full_worker.replace(
    ".rc_full_gem_step2_reaction_batch_worker", ".rc_step2_reaction_batch_worker"
).replace(
    "Malformed full-GEM Step 2 reaction-batch task.", "Malformed Step 2 reaction-batch task."
).replace(
    "regcompass_full_gem_step2_compact_payload_v1", "regcompass_step2_compact_payload_v1"
).replace(
    "Malformed full-GEM Step 2 compact payload.", "Malformed Step 2 compact payload."
).replace(
    "Full-GEM Step 2 reaction-batch rows are absent from the payload.",
    "Step 2 reaction-batch rows are absent from the compact payload."
).replace(
    "Full-GEM Step 2 compact payload lacks the required LP model state.",
    "Step 2 compact payload lacks the required LP model state."
).replace(
    "Full-GEM Step 2 payload reaction order differs from its model.",
    "Step 2 compact payload reaction order differs from its union GEM."
).replace(
    "Full-GEM Step 2 compact payload penalties are not aligned.",
    "Step 2 compact payload penalties are not aligned."
).replace(
    "Full-GEM paired RNA-control penalties are not aligned.",
    "Paired RNA-control penalties are not aligned."
).replace(
    "A full-GEM target reaction is absent from its shared model.",
    "A target reaction is absent from its cell-type union GEM."
).replace(
    "module_id = rep(NA_character_, length(units)),",
    "module_id = rep(\"CELLTYPE_MEDIUM_UNION_GEM\", length(units)),\n      cell_type = rep(entry$cell_type, length(units)),"
).replace(
    "condition = rep(entry$condition, length(units)),",
    "condition = rep(\"all\", length(units)),"
).replace(
    'primary$feasible, "ok", "medium_directionally_infeasible"',
    'primary$feasible, "ok", "structurally_infeasible"'
).replace(
    "vmax_reused_from_shared_cache = rep(TRUE, length(units)),",
    "vmax_reused_from_celltype_cache = rep(TRUE, length(units)),"
).replace(
    '"directional_reaction_x_all_metacells", length(units)',
    '"directional_reaction_x_matching_metacells", length(units)'
)
# Cell-type entries do not have condition; control diagnostics are fine without it.
text = replace_regex(
    text,
    r'\.rc_step2_reaction_batch_worker <- function\(task\) \{.*?\n}\n\n(?=\.rc_run_celltype_microcompass_engine_reaction_core)',
    cell_worker,
    "replace celltype worker"
)
text = replace_once(
    text,
    "    BPPARAM = NULL,\n    model_cache_override = NULL) {\n  unit <- match.arg(unit)",
    "    BPPARAM = NULL,\n    model_cache_override = NULL,\n    control_layer1 = NULL) {\n  unit <- match.arg(unit)",
    "celltype core control arg"
)
text = replace_once(
    text,
    '''  matrices <- rc_layer2_unit_matrices(\n    layer1,\n    if (identical(unit, "metacell")) "metacell" else "sample_celltype",\n    sample_col, celltype_col, condition_col\n  )\n  unit_meta <- matrices$unit_meta''',
    '''  matrices <- rc_layer2_unit_matrices(\n    layer1,\n    if (identical(unit, "metacell")) "metacell" else "sample_celltype",\n    sample_col, celltype_col, condition_col\n  )\n  control_matrices <- if (!is.null(control_layer1)) {\n    value <- rc_layer2_unit_matrices(\n      control_layer1,\n      if (identical(unit, "metacell")) "metacell" else "sample_celltype",\n      sample_col, celltype_col, condition_col\n    )\n    if (!identical(dimnames(value$reaction_expression),\n                   dimnames(matrices$reaction_expression)) ||\n        !identical(value$unit_meta, matrices$unit_meta)) {\n      stop("Paired RNA-control units do not align with primary Layer 2 units.",\n           call. = FALSE)\n    }\n    value\n  } else NULL\n  unit_meta <- matrices$unit_meta''',
    "celltype control matrices"
)
text = replace_once(
    text,
    '''  penalties <- rc_compute_multiome_penalty(\n    rc_align_reaction_expression(\n      matrices$reaction_expression, all_reactions, NA_real_\n    ),\n    reaction_roles = gem$reaction_roles\n  )\n  vmax_cache <- .rc_build_microcompass_vmax_cache(''',
    '''  penalties <- rc_compute_multiome_penalty(\n    rc_align_reaction_expression(\n      matrices$reaction_expression, all_reactions, NA_real_\n    ),\n    reaction_roles = gem$reaction_roles\n  )\n  control_penalties <- if (!is.null(control_matrices)) {\n    rc_compute_multiome_penalty(\n      rc_align_reaction_expression(\n        control_matrices$reaction_expression, all_reactions, NA_real_\n      ),\n      reaction_roles = gem$reaction_roles\n    )\n  } else NULL\n  vmax_cache <- .rc_build_microcompass_vmax_cache(''',
    "celltype control penalties"
)
text = replace_once(
    text,
    '''  feasible <- evaluated <- matrix(\n    FALSE, length(row_ids), length(units),\n    dimnames = list(row_ids, units)\n  )\n\n  checkpoint_root''',
    '''  feasible <- evaluated <- matrix(\n    FALSE, length(row_ids), length(units),\n    dimnames = list(row_ids, units)\n  )\n  control_penalty <- if (!is.null(control_penalties)) penalty else NULL\n  control_feasible <- if (!is.null(control_penalties)) feasible else NULL\n  control_evaluated <- if (!is.null(control_penalties)) evaluated else NULL\n\n  checkpoint_root''',
    "celltype control output matrices"
)
text = replace_once(
    text,
    '''      flux_threshold = flux_threshold,\n      payload_dir = payload_dir\n    )''',
    '''      flux_threshold = flux_threshold,\n      payload_dir = payload_dir,\n      control_penalties = control_penalties\n    )''',
    "celltype payload control pass"
)
text = replace_once(
    text,
    '''  penalties$penalty <- NULL\n  rm(vmax_cache, matrices, all_reactions)''',
    '''  penalties$penalty <- NULL\n  if (!is.null(control_penalties)) control_penalties$penalty <- NULL\n  rm(vmax_cache, matrices, control_matrices, all_reactions)''',
    "celltype cleanup control"
)
text = replace_once(
    text,
    '''  diagnostics <- vector("list", length(checkpoint_files))\n  step2_engine_metrics <- vector("list", length(checkpoint_files))\n  observed_rows <- character(length(checkpoint_files))''',
    '''  diagnostics <- vector("list", length(checkpoint_files))\n  step2_engine_metrics <- vector("list", length(checkpoint_files))\n  control_diagnostics <- vector("list", length(checkpoint_files))\n  control_step2_engine_metrics <- vector("list", length(checkpoint_files))\n  control_reused <- logical(length(checkpoint_files))\n  observed_rows <- character(length(checkpoint_files))''',
    "celltype control collector init"
)
text = replace_once(
    text,
    '''    step2_engine_metrics[[i]] <- data.frame(\n      row_id = row_id,\n      engine = as.character(metrics$engine %||% "unknown"),\n      n_solves = as.integer(metrics$n_solves %||% 0L),\n      n_objective_updates = as.integer(\n        metrics$n_objective_updates %||% 0L\n      ),\n      n_fallback = as.integer(metrics$n_fallback %||% 0L),\n      stringsAsFactors = FALSE\n    )\n    rm(result)''',
    '''    step2_engine_metrics[[i]] <- data.frame(\n      row_id = row_id,\n      engine = as.character(metrics$engine %||% "unknown"),\n      n_solves = as.integer(metrics$n_solves %||% 0L),\n      n_objective_updates = as.integer(\n        metrics$n_objective_updates %||% 0L\n      ),\n      n_fallback = as.integer(metrics$n_fallback %||% 0L),\n      stringsAsFactors = FALSE\n    )\n    if (!is.null(control_penalty)) {\n      ctrl <- result$control\n      if (!is.list(ctrl) ||\n          !identical(names(ctrl$penalty), as.character(result$units))) {\n        stop("A paired RNA-control checkpoint is malformed.", call. = FALSE)\n      }\n      control_penalty[row_id, result$units] <- ctrl$penalty\n      control_feasible[row_id, result$units] <- ctrl$feasible\n      control_evaluated[row_id, result$units] <- ctrl$evaluated\n      control_diagnostics[[i]] <- ctrl$diagnostics\n      control_reused[[i]] <- isTRUE(ctrl$reused_from_primary)\n      control_metrics <- ctrl$engine_metrics %||% list()\n      control_step2_engine_metrics[[i]] <- data.frame(\n        row_id = row_id,\n        engine = as.character(control_metrics$engine %||% "not_run"),\n        n_solves = as.integer(control_metrics$n_solves %||% 0L),\n        n_objective_updates = as.integer(\n          control_metrics$n_objective_updates %||% 0L\n        ),\n        n_fallback = as.integer(control_metrics$n_fallback %||% 0L),\n        reused_from_primary = isTRUE(ctrl$reused_from_primary),\n        stringsAsFactors = FALSE\n      )\n    }\n    rm(result)''',
    "celltype collect control"
)
text = replace_once(
    text,
    '''  score <- rc_compass_score_from_penalty(penalty, feasible)\n  directions <- unique(do.call(rbind, lapply(model_cache, function(entry) {''',
    '''  score <- rc_compass_score_from_penalty(penalty, feasible)\n  score_rna_only <- if (!is.null(control_penalty)) {\n    rc_compass_score_from_penalty(control_penalty, control_feasible)\n  } else NULL\n  directions <- unique(do.call(rbind, lapply(model_cache, function(entry) {''',
    "celltype control score"
)
text = replace_once(
    text,
    '''    evaluated = evaluated,\n    target_direction = directions,''',
    '''    evaluated = evaluated,\n    penalty_rna_only = control_penalty,\n    score_rna_only = score_rna_only,\n    feasible_rna_only = control_feasible,\n    evaluated_rna_only = control_evaluated,\n    lp_diagnostics_rna_only = .rc_bind_frames_fill(control_diagnostics),\n    step2_engine_metrics_rna_only =\n      .rc_bind_frames_fill(control_step2_engine_metrics),\n    rna_control_model_identical_reuse = control_reused,\n    target_direction = directions,''',
    "celltype return control"
)
text = replace_once(
    text,
    '''      step2_solver_reuse = paste(\n        "one persistent HiGHS model per directional reaction reused across",\n        "all matching metacells; one-shot fallback for unsupported backends"\n      ),''',
    '''      step2_solver_reuse = paste(\n        "one prepared target template with independent persistent HiGHS",\n        "streams for primary and RNA-only objectives across matching metacells"\n      ),\n      paired_primary_rna_control = !is.null(control_penalty),\n      rna_control_vmax_solve_count = if (!is.null(control_penalty)) 0L else NA_integer_,\n      rna_control_identical_model_reuse = sum(control_reused),''',
    "celltype params control"
)
write(path, text)

# Stage-level wiring: create RNA-only Layer 1 once, pass it into the primary engine,
# remove the second complete engine invocation, and reconstruct the compatibility path.
path = "R/step_layer2.R"
text = read(path)
text = replace_once(
    text,
    '''    "parallel", "workers", "penalty_weights"\n  ))''',
    '''    "parallel", "workers", "penalty_weights", "control_layer1"\n  ))''',
    "reserve control_layer1"
)
text = replace_once(
    text,
    '''  defaults <- list(\n    layer1 = layer1,''',
    '''  rna_expression <- layer1$reaction_expression_rna_only\n  if (!is.numeric(rna_expression) || is.null(dim(rna_expression)) ||\n      !identical(dimnames(rna_expression),\n                 dimnames(layer1$reaction_expression))) {\n    stop("Layer 1 RNA-only control is missing or misaligned.", call. = FALSE)\n  }\n  control_layer1 <- layer1\n  control_layer1$reaction_expression <- rna_expression\n\n  defaults <- list(\n    layer1 = layer1,\n    control_layer1 = control_layer1,''',
    "stage create control layer"
)
# Replace old post-primary full control run block.
old_pattern = r'''  control_model_cache <- answer\$shared_model_cache.*?  answer\$comparison_paths <- list\(rna_only = rna_only\)'''
new_block = r'''  if (is.null(answer$penalty_rna_only) || is.null(answer$score_rna_only) ||
      is.null(answer$feasible_rna_only) || is.null(answer$evaluated_rna_only)) {
    stop("Paired Layer 2 did not return the RNA-only control matrices.",
         call. = FALSE)
  }
  rna_only <- answer
  rna_only$penalty <- answer$penalty_rna_only
  rna_only$score <- answer$score_rna_only
  rna_only$feasible <- answer$feasible_rna_only
  rna_only$evaluated <- answer$evaluated_rna_only
  rna_only$lp_diagnostics <- answer$lp_diagnostics_rna_only
  rna_only$step2_engine_metrics <- answer$step2_engine_metrics_rna_only
  rna_only$shared_model_cache <- NULL
  .rc_assert_layer2_shared_contract(answer, rna_only, "RNA-only control")

  answer$schema_version <- "regcompass_regulatory_layer2_v3"
  answer$comparison_contract <- list(
    primary = "penalty",
    rna_control = "penalty_rna_only",
    nonestimable_edge_policy =
      "coefficient_NA_and_zero_realized_penalty_contribution",
    exact_shared_structure = TRUE,
    structural_model_contract = answer$structural_model_contract,
    directional_vmax_contract =
      "computed_once_and_shared_by_paired_primary_rna_control",
    rna_control_vmax_solve_count = 0L,
    paired_step2_dispatch = TRUE,
    independent_solver_streams = TRUE,
    exact_identical_model_objective_reuse = TRUE,
    effect_size_basis = "penalty / (omega * vmax)",
    ecdf_effect_size_eligible = FALSE
  )
  answer$comparison_paths <- list(rna_only = rna_only)'''
text = replace_regex(text, old_pattern, new_block, "replace second Layer2 run")
# Progress text now describes paired execution; keep legacy control wrapper unused for compatibility.
text = replace_once(
    text,
    '    "starting structural construction and primary multiome scoring"',
    '    "starting structural construction and paired primary/RNA-only scoring"',
    "paired progress start"
)
text = replace_once(
    text,
    '''  .rc_layer2_overall_event(\n    "primary_scoring_complete", 5L,\n    "primary multiome Step 2 scoring completed"\n  )''',
    '''  .rc_layer2_overall_event(\n    "primary_scoring_complete", 5L,\n    "paired primary multiome and RNA-only Step 2 scoring completed"\n  )\n  if (!is.null(args$control_layer1)) {\n    .rc_layer2_overall_event(\n      "rna_control_complete", 8L,\n      "RNA-only Step 2 completed in the same target dispatch with shared Vmax"\n    )\n  }''',
    "paired progress complete"
)
write(path, text)

# Standalone numerical equivalence regression.
test = r'''options(warn = 2)
source("R/00_utils.R", local = FALSE)
source("R/lp_solver.R", local = FALSE)
source("R/layer2_parallel_runtime.R", local = FALSE)
source("R/microcompass_vmax_cache.R", local = FALSE)
source("R/microcompass_engine.R", local = FALSE)
source("R/celltype_microcompass_reaction_parallel.R", local = FALSE)

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
  stopifnot(
    isTRUE(all.equal(paired$penalty, primary_legacy$penalty,
                     tolerance = 1e-12)),
    isTRUE(all.equal(paired$feasible, primary_legacy$feasible)),
    isTRUE(all.equal(paired$evaluated, primary_legacy$evaluated)),
    isTRUE(all.equal(paired$control$penalty, control_legacy$penalty,
                     tolerance = 1e-12)),
    isTRUE(all.equal(paired$control$feasible, control_legacy$feasible)),
    isTRUE(all.equal(paired$control$evaluated, control_legacy$evaluated)),
    paired$engine_metrics$n_solves == length(units),
    paired$control$engine_metrics$n_solves == length(units),
    !isTRUE(paired$control$reused_from_primary)
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
}

check_worker(.rc_full_gem_step2_reaction_batch_worker, TRUE)
check_worker(.rc_step2_reaction_batch_worker, FALSE)

stage <- paste(readLines("R/step_layer2.R", warn = FALSE), collapse = "\n")
stopifnot(
  grepl("control_layer1 = control_layer1", stage, fixed = TRUE),
  grepl("paired_step2_dispatch = TRUE", stage, fixed = TRUE),
  !grepl("rna_only <- run_control(", stage, fixed = TRUE)
)

cat("paired Layer 2 primary/RNA-control numerical equivalence checks passed\n")
'''
write("tests/layer2-paired-control-check.R", test)

# testthat wrapper is intentionally structural; the standalone script carries numeric HiGHS checks.
testthat = r'''test_that("Layer 2 primary and RNA control share one target dispatch", {
  stage <- paste(readLines(file.path("R", "step_layer2.R"), warn = FALSE),
                 collapse = "\n")
  full <- paste(readLines(file.path("R", "microcompass_engine.R"), warn = FALSE),
                collapse = "\n")
  cell <- paste(readLines(
    file.path("R", "celltype_microcompass_reaction_parallel.R"),
    warn = FALSE
  ), collapse = "\n")
  expect_match(stage, "control_layer1 = control_layer1", fixed = TRUE)
  expect_match(stage, "paired_step2_dispatch = TRUE", fixed = TRUE)
  expect_false(grepl("rna_only <- run_control(", stage, fixed = TRUE))
  expect_match(full, "control_penalties = control_penalties", fixed = TRUE)
  expect_match(cell, "control_penalties = control_penalties", fixed = TRUE)
  expect_match(full, ".rc_compass_step2_route_solve", fixed = TRUE)
  expect_match(cell, ".rc_compass_step2_route_solve", fixed = TRUE)
})
'''
write("tests/testthat/test-layer2-paired-control.R", testthat)

print("Layer 2 paired-control patch applied")
