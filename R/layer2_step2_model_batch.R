# Exact model-batch acceleration for Layer 2 COMPASS Step 2.
#
# Cross-target solver reuse is restricted to COMPASS Step 2. Production Step 2
# consumes the scalar optimum/status, not a particular primal support pattern,
# so warm-start/basis reuse cannot change the biological estimand. CORDA2 keeps
# its target-local solver ownership because CORDA2 does consume primal support
# for structural reconstruction.

.rc_step2_batch_engine_metrics <- function(engine) {
  base <- .rc_compass_step2_engine_metrics(engine)
  c(
    base,
    list(
      n_bound_updates = as.integer(
        if (is.list(engine)) engine$batch_n_bound_updates %||% 0L else 0L
      ),
      n_target_switches = as.integer(
        if (is.list(engine)) engine$batch_n_target_switches %||% 0L else 0L
      )
    )
  )
}

.rc_step2_batch_metric_delta <- function(after, before) {
  fields <- c(
    "n_solves", "n_objective_updates", "n_fallback",
    "n_bound_updates", "n_target_switches"
  )
  out <- lapply(fields, function(field) {
    max(
      0L,
      as.integer(after[[field]] %||% 0L) -
        as.integer(before[[field]] %||% 0L)
    )
  })
  names(out) <- fields
  c(list(engine = as.character(after$engine %||% "not_run")), out)
}

.rc_step2_batch_accumulate_metric <- function(total, delta) {
  fields <- c(
    "n_solves", "n_objective_updates", "n_fallback",
    "n_bound_updates", "n_target_switches"
  )
  if (is.null(total)) {
    total <- c(
      list(engine = as.character(delta$engine %||% "not_run")),
      stats::setNames(as.list(rep.int(0L, length(fields))), fields),
      list(
        n_objective_change_events = 0L,
        n_bound_change_events = 0L
      )
    )
  }
  if (!identical(as.character(delta$engine %||% "not_run"), "not_run")) {
    total$engine <- as.character(delta$engine)
  }
  for (field in fields) {
    total[[field]] <- as.integer(total[[field]] %||% 0L) +
      as.integer(delta[[field]] %||% 0L)
  }
  if (as.integer(delta$n_objective_updates %||% 0L) > 0L) {
    total$n_objective_change_events <-
      as.integer(total$n_objective_change_events %||% 0L) + 1L
  }
  if (as.integer(delta$n_bound_updates %||% 0L) > 0L) {
    total$n_bound_change_events <-
      as.integer(total$n_bound_change_events %||% 0L) + 1L
  }
  total
}

.rc_step2_batch_early_result <- function(vmax_result) {
  list(
    feasible = FALSE,
    penalty = NA_real_,
    vmax = as.numeric(vmax_result$vmax),
    solver_status = as.character(vmax_result$status),
    step1_status = as.character(vmax_result$status),
    step2_status = "not_run",
    solver_backend = "not_run",
    flux = numeric()
  )
}

.rc_step2_batch_compile_shared <- function(payload, row_ids) {
  model <- payload$model
  reactions <- as.character(payload$reactions)
  row_ids <- as.character(row_ids)
  feasible_rows <- row_ids[vapply(
    payload$vmax[row_ids],
    function(value) isTRUE(value$feasible),
    logical(1)
  )]
  if (!length(feasible_rows)) return(NULL)

  first_row <- feasible_rows[[1L]]
  first_entry <- payload$entries[[first_row]]
  first <- .rc_compass_step2_prepare(
    S = model$S,
    lb = model$lb,
    ub = model$ub,
    target_reaction = first_entry$reaction_id,
    vmax_result = payload$vmax[[first_row]],
    target_direction = first_entry$target_direction,
    omega = payload$omega,
    flux_threshold = payload$flux_threshold
  )
  if (!isTRUE(first$runnable)) {
    stop("A feasible directional Vmax did not produce a runnable Step 2 model.",
         call. = FALSE)
  }

  template <- first$template
  signed_lb <- rc_align_bound(
    model$lb, reactions, default = -1000, name = "lb"
  )
  signed_ub <- rc_align_bound(
    model$ub, reactions, default = 1000, name = "ub"
  )
  reaction_index <- template$direction_reaction_index
  direction_sign <- template$direction_sign
  base_lower <- ifelse(
    direction_sign > 0,
    pmax(signed_lb[reaction_index], 0),
    pmax(-signed_ub[reaction_index], 0)
  )
  base_upper <- ifelse(
    direction_sign > 0,
    signed_ub[reaction_index],
    -signed_lb[reaction_index]
  )
  base_lower <- as.numeric(base_lower)
  base_upper <- as.numeric(base_upper)
  if (length(base_lower) != template$n_variables ||
      length(base_upper) != template$n_variables ||
      any(!is.finite(base_lower)) || any(!is.finite(base_upper)) ||
      any(base_lower > base_upper)) {
    stop("The shared directional Step 2 base bounds are malformed.",
         call. = FALSE)
  }

  # S_dir and directional variable identity depend only on the completed GEM.
  # Reset the first target-specific template to the exact target-independent
  # completed-GEM bounds before one persistent solver is created.
  template$lb <- base_lower
  template$ub <- base_upper
  template$vmax <- NA_real_
  template$step1_status <- "shared_model_template"
  template$target_reaction <- NA_character_
  template$target_direction <- NA_character_
  template$target_variable_index <- integer()
  template$opposite_variable_index <- integer()
  template$required_flux <- NA_real_
  template$formulation <-
    "compass_directional_nonnegative_exact_model_batch_v1"

  list(
    template = template,
    base_lower = base_lower,
    base_upper = base_upper
  )
}

.rc_step2_batch_target_spec <- function(shared, payload, row_id) {
  reactions <- as.character(payload$reactions)
  entry <- payload$entries[[row_id]]
  vmax_result <- payload$vmax[[row_id]]
  reaction <- as.character(entry$reaction_id)
  direction <- as.character(entry$target_direction)

  if (length(reaction) != 1L || is.na(reaction) || !nzchar(reaction) ||
      !reaction %in% reactions) {
    stop("A Step 2 batch target reaction is absent from its shared model.",
         call. = FALSE)
  }
  if (length(direction) != 1L ||
      !direction %in% c("forward", "reverse")) {
    stop("A Step 2 batch target direction must be forward or reverse.",
         call. = FALSE)
  }
  if (!is.list(vmax_result) ||
      !all(c("feasible", "vmax", "status") %in% names(vmax_result))) {
    stop("A Step 2 batch target lacks a directional Vmax result.",
         call. = FALSE)
  }

  target_index <- match(reaction, reactions)
  if (!isTRUE(vmax_result$feasible)) {
    return(list(
      runnable = FALSE,
      target_index = as.integer(target_index),
      result = .rc_step2_batch_early_result(vmax_result)
    ))
  }
  if (is.null(shared)) {
    stop("A feasible Step 2 target requires a shared directional model.",
         call. = FALSE)
  }

  vmax <- as.numeric(vmax_result$vmax)
  if (length(vmax) != 1L || !is.finite(vmax) ||
      vmax < payload$flux_threshold) {
    stop("Cached directional vmax is not a positive feasible value.",
         call. = FALSE)
  }
  target_sign <- if (identical(direction, "forward")) 1 else -1
  target_variable <- which(
    shared$template$direction_reaction_index == target_index &
      shared$template$direction_sign == target_sign
  )
  if (length(target_variable) != 1L) {
    stop(
      "Cached directional vmax is feasible but the requested target direction ",
      "is absent from the shared Step 2 directional model.",
      call. = FALSE
    )
  }
  opposite_variable <- which(
    shared$template$direction_reaction_index == target_index &
      shared$template$direction_sign == -target_sign
  )
  required_flux <- as.numeric(payload$omega) * vmax
  bound_tolerance <- max(
    payload$flux_threshold,
    32 * .Machine$double.eps * max(
      1, abs(required_flux), abs(shared$base_upper[[target_variable]])
    )
  )
  if (required_flux >
      shared$base_upper[[target_variable]] + bound_tolerance) {
    stop(
      "Cached directional vmax is inconsistent with the shared completed-GEM ",
      "target bound.", call. = FALSE
    )
  }
  target_lower <- max(
    shared$base_lower[[target_variable]],
    min(required_flux, shared$base_upper[[target_variable]])
  )
  if (length(opposite_variable) &&
      shared$base_lower[[opposite_variable]] > bound_tolerance) {
    stop(
      "The opposite target direction has a positive compulsory lower bound; ",
      "this contradicts the feasible signed directional Vmax contract.",
      call. = FALSE
    )
  }

  list(
    runnable = TRUE,
    target_index = as.integer(target_index),
    target_variable_index = as.integer(target_variable),
    opposite_variable_index = as.integer(opposite_variable),
    target_lower = as.numeric(target_lower),
    required_flux = as.numeric(required_flux),
    vmax = vmax,
    step1_status = as.character(vmax_result$status),
    target_reaction = reaction,
    target_direction = direction,
    bound_indices = as.integer(unique(c(target_variable, opposite_variable)))
  )
}

.rc_step2_batch_new_engine <- function(shared, solver) {
  if (is.null(shared)) return(NULL)
  engine <- .rc_compass_step2_new_engine(
    shared$template,
    solver,
    persistent_required = identical(as.character(solver), "highs")
  )
  if (identical(as.character(solver), "highs") &&
      !"hi_solver_set_variable_bounds" %in% getNamespaceExports("highs")) {
    .rc_compass_step2_release_engine(engine)
    stop(
      "Model-batch Layer 2 requires HiGHS sparse variable-bound updates. ",
      "Install a current `highs` package (>= 1.12.0-1).",
      call. = FALSE
    )
  }
  engine$batch_base_lower <- shared$base_lower
  engine$batch_base_upper <- shared$base_upper
  engine$batch_current_target_indices <- integer()
  engine$batch_n_bound_updates <- 0L
  engine$batch_n_target_switches <- 0L
  engine
}

.rc_step2_batch_apply_target <- function(engine, spec) {
  if (!isTRUE(spec$runnable)) return(engine)
  if (!is.list(engine) || is.null(engine$template)) {
    stop("A runnable Step 2 target requires an initialized batch engine.",
         call. = FALSE)
  }

  # Include both the previous and requested target indices. Starting candidate
  # bounds from immutable completed-GEM base bounds guarantees that no target-
  # specific lower/upper bound can leak into the next LP.
  previous <- as.integer(engine$batch_current_target_indices %||% integer())
  requested <- as.integer(spec$bound_indices)
  candidate <- sort(unique(c(previous, requested)))
  lower <- engine$batch_base_lower[candidate]
  upper <- engine$batch_base_upper[candidate]

  target_position <- match(spec$target_variable_index, candidate)
  lower[[target_position]] <- spec$target_lower
  if (length(spec$opposite_variable_index)) {
    opposite_position <- match(spec$opposite_variable_index, candidate)
    lower[[opposite_position]] <- 0
    upper[[opposite_position]] <- 0
  }
  changed <- candidate[
    engine$template$lb[candidate] != lower |
      engine$template$ub[candidate] != upper
  ]

  persistent_error <- NULL
  if (length(changed) &&
      identical(engine$type, "highs_persistent_cpp") &&
      !isTRUE(engine$persistent_disabled)) {
    positions <- match(changed, candidate)
    persistent_error <- tryCatch({
      .rc_microcompass_highs_call(
        "hi_solver_set_variable_bounds", engine$pointer,
        index = as.integer(changed - 1L),
        lower = lower[positions],
        upper = upper[positions]
      )
      NULL
    }, error = function(error) error)
  }
  if (inherits(persistent_error, "error")) {
    engine$persistent_message <- conditionMessage(persistent_error)
    if (isTRUE(engine$persistent_required)) {
      .rc_compass_step2_release_engine(engine)
      stop(
        "Persistent HiGHS target-bound update failed inside model-batch ",
        "Layer 2: ", engine$persistent_message,
        call. = FALSE
      )
    }
    engine <- .rc_compass_step2_release_engine(engine)
    engine$type <- "one_shot"
    engine$persistent_disabled <- TRUE
    engine$n_fallback <- as.integer(engine$n_fallback %||% 0L) + 1L
  }

  engine$template$lb[candidate] <- lower
  engine$template$ub[candidate] <- upper
  engine$template$vmax <- spec$vmax
  engine$template$step1_status <- spec$step1_status
  engine$template$target_reaction <- spec$target_reaction
  engine$template$target_direction <- spec$target_direction
  engine$template$target_variable_index <- spec$target_variable_index
  engine$template$opposite_variable_index <- spec$opposite_variable_index
  engine$template$required_flux <- spec$required_flux
  engine$batch_current_target_indices <- requested
  engine$batch_n_bound_updates <-
    as.integer(engine$batch_n_bound_updates %||% 0L) + length(changed)
  engine$batch_n_target_switches <-
    as.integer(engine$batch_n_target_switches %||% 0L) + 1L
  engine
}

.rc_step2_batch_route_result_template <- function(
    spec, penalty_matrix, units) {
  n_units <- length(units)
  target_available <- is.finite(
    penalty_matrix[spec$target_index, , drop = TRUE]
  )
  names(target_available) <- units
  list(
    penalty = stats::setNames(rep(NA_real_, n_units), units),
    vmax = stats::setNames(rep(NA_real_, n_units), units),
    feasible = stats::setNames(rep(FALSE, n_units), units),
    evaluated = stats::setNames(rep(FALSE, n_units), units),
    solver_status = stats::setNames(rep(NA_character_, n_units), units),
    solver_backend = stats::setNames(rep(NA_character_, n_units), units),
    step1_status = stats::setNames(rep(NA_character_, n_units), units),
    step2_status = stats::setNames(rep(NA_character_, n_units), units),
    target_available = target_available
  )
}

.rc_step2_batch_route_solve <- function(
    specs, engine, penalty_matrix, evidence, units,
    reuse_mask = NULL, reuse_results = NULL) {
  penalty_matrix <- as.matrix(penalty_matrix)
  units <- as.character(units)
  row_ids <- names(specs)
  n_units <- length(units)
  if (is.null(row_ids) || any(!nzchar(row_ids)) ||
      !identical(colnames(penalty_matrix), units)) {
    stop("Model-batch Step 2 route inputs are not aligned.", call. = FALSE)
  }
  if (length(evidence$all_finite) != n_units ||
      length(evidence$fraction) != n_units ||
      length(evidence$unavailable) != n_units) {
    stop("Model-batch Step 2 penalty evidence summary is malformed.",
         call. = FALSE)
  }
  reuse_mask <- if (is.null(reuse_mask)) {
    rep(FALSE, n_units)
  } else {
    as.logical(reuse_mask)
  }
  if (length(reuse_mask) != n_units || anyNA(reuse_mask)) {
    stop("Model-batch Step 2 reuse mask is not aligned to units.",
         call. = FALSE)
  }
  if (any(reuse_mask) &&
      (is.null(reuse_results) || !all(row_ids %in% names(reuse_results)))) {
    stop("Model-batch Step 2 exact-reuse source is malformed.",
         call. = FALSE)
  }

  results <- lapply(specs, function(spec) {
    .rc_step2_batch_route_result_template(spec, penalty_matrix, units)
  })
  metrics <- lapply(specs, function(spec) {
    list(
      engine = if (isTRUE(spec$runnable) && is.list(engine)) {
        as.character(engine$type %||% "one_shot")
      } else {
        "not_run"
      },
      n_solves = 0L,
      n_objective_updates = 0L,
      n_fallback = 0L,
      n_bound_updates = 0L,
      n_target_switches = 0L,
      n_objective_change_events = 0L,
      n_bound_change_events = 0L
    )
  })
  names(results) <- names(metrics) <- row_ids
  route_objective_events <- 0L
  route_bound_events <- 0L

  # Unit-major traversal is exact: the objective is unit-specific, while the
  # target changes only variable bounds. One unit objective can therefore be
  # reused across every target in this structural-model batch.
  for (i in seq_along(units)) {
    if (isTRUE(reuse_mask[[i]])) {
      for (row_id in row_ids) {
        source <- reuse_results[[row_id]]
        required <- c(
          "penalty", "vmax", "feasible", "evaluated", "solver_status",
          "step1_status", "step2_status", "target_available"
        )
        if (!is.list(source) || !all(required %in% names(source))) {
          stop("A model-batch Step 2 reuse result is malformed.",
               call. = FALSE)
        }
        results[[row_id]]$penalty[[i]] <- source$penalty[[i]]
        results[[row_id]]$vmax[[i]] <- source$vmax[[i]]
        results[[row_id]]$feasible[[i]] <- source$feasible[[i]]
        results[[row_id]]$evaluated[[i]] <- source$evaluated[[i]]
        results[[row_id]]$solver_status[[i]] <- source$solver_status[[i]]
        results[[row_id]]$solver_backend[[i]] <-
          "reused_identical_primary_objective"
        results[[row_id]]$step1_status[[i]] <- source$step1_status[[i]]
        results[[row_id]]$step2_status[[i]] <- source$step2_status[[i]]
      }
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

    for (row_id in row_ids) {
      spec <- specs[[row_id]]
      before <- .rc_step2_batch_engine_metrics(engine)
      if (isTRUE(spec$runnable)) {
        engine <- .rc_step2_batch_apply_target(engine, spec)
        solved <- .rc_compass_step2_engine_solve(
          engine,
          solver_penalty,
          return_solution = FALSE,
          trusted_aligned = TRUE
        )
        engine <- solved$engine
        fit <- .rc_compass_step2_result(
          engine$template,
          solved$answer,
          require_solution = FALSE
        )
      } else {
        fit <- spec$result
        solved <- NULL
      }
      after <- .rc_step2_batch_engine_metrics(engine)
      delta <- .rc_step2_batch_metric_delta(after, before)
      metrics[[row_id]] <-
        .rc_step2_batch_accumulate_metric(metrics[[row_id]], delta)
      if (delta$n_objective_updates > 0L) {
        route_objective_events <- route_objective_events + 1L
      }
      if (delta$n_bound_updates > 0L) {
        route_bound_events <- route_bound_events + 1L
      }

      available <- results[[row_id]]$target_available[[i]]
      results[[row_id]]$penalty[[i]] <-
        if (available) fit$penalty else NA_real_
      results[[row_id]]$vmax[[i]] <- fit$vmax
      results[[row_id]]$feasible[[i]] <- isTRUE(fit$feasible)
      results[[row_id]]$evaluated[[i]] <-
        isTRUE(fit$feasible) && available
      results[[row_id]]$solver_status[[i]] <-
        as.character(fit$solver_status)
      results[[row_id]]$solver_backend[[i]] <-
        as.character(fit$solver_backend %||% "unknown")
      results[[row_id]]$step1_status[[i]] <-
        as.character(fit$step1_status)
      results[[row_id]]$step2_status[[i]] <-
        as.character(fit$step2_status)
      rm(fit, solved, before, after, delta)
    }
    rm(unit_penalty, solver_penalty)
  }

  list(
    engine = engine,
    results = results,
    metrics = metrics,
    batch_metrics = list(
      n_solves = as.integer(sum(vapply(
        metrics, function(x) x$n_solves, integer(1)
      ))),
      n_objective_change_events = as.integer(route_objective_events),
      n_bound_change_events = as.integer(route_bound_events),
      n_target_switches = as.integer(sum(vapply(
        metrics, function(x) x$n_target_switches, integer(1)
      ))),
      traversal = "unit_then_directional_target",
      shared_model_solver = TRUE
    )
  )
}

.rc_step2_batch_validate_worker_payload <- function(task, full_gem = FALSE) {
  if (!is.list(task) ||
      !all(c("payload_file", "row_ids", "checkpoint_dir") %in% names(task))) {
    stop("Malformed Step 2 model-batch task.", call. = FALSE)
  }
  payload <- readRDS(task$payload_file)
  expected_schema <- if (isTRUE(full_gem)) {
    "regcompass_full_gem_step2_compact_payload_v1"
  } else {
    "regcompass_step2_compact_payload_v1"
  }
  if (!is.list(payload) ||
      !identical(payload$schema_version, expected_schema)) {
    stop("Malformed Step 2 model-batch compact payload.", call. = FALSE)
  }
  row_ids <- as.character(task$row_ids)
  if (!length(row_ids) ||
      !all(row_ids %in% names(payload$entries)) ||
      !all(row_ids %in% names(payload$vmax))) {
    stop("Step 2 model-batch target rows are absent from the payload.",
         call. = FALSE)
  }
  model <- payload$model
  if (!is.list(model) || is.null(model$S) || is.null(model$lb) ||
      is.null(model$ub)) {
    stop("Step 2 model-batch payload lacks the required LP model state.",
         call. = FALSE)
  }
  if (!identical(colnames(model$S), as.character(payload$reactions))) {
    stop("Step 2 model-batch reaction order differs from its shared model.",
         call. = FALSE)
  }
  if (!identical(colnames(payload$penalty), as.character(payload$units)) ||
      !identical(rownames(payload$penalty), as.character(payload$reactions))) {
    stop("Step 2 model-batch penalties are not aligned.", call. = FALSE)
  }
  paired_control <- !is.null(payload$control_penalty)
  if (paired_control &&
      !identical(dimnames(payload$control_penalty), dimnames(payload$penalty))) {
    stop("Paired RNA-control penalties are not aligned.", call. = FALSE)
  }
  list(
    payload = payload,
    row_ids = row_ids,
    model = model,
    paired_control = paired_control
  )
}

.rc_step2_batch_run_routes <- function(payload, row_ids, paired_control) {
  units <- as.character(payload$units)
  evidence <- payload$penalty_evidence %||%
    .rc_step2_penalty_evidence_stats(payload$penalty)
  control_evidence <- if (paired_control) {
    payload$control_penalty_evidence %||%
      .rc_step2_penalty_evidence_stats(payload$control_penalty)
  } else {
    NULL
  }
  reuse_mask <- if (paired_control) {
    value <- as.logical(payload$control_identical)
    if (length(value) != length(units) || anyNA(value)) {
      stop("Paired RNA-control exact-reuse mask is not aligned to units.",
           call. = FALSE)
    }
    names(value) <- units
    value
  } else {
    stats::setNames(rep(FALSE, length(units)), units)
  }

  shared <- .rc_step2_batch_compile_shared(payload, row_ids)
  specs <- lapply(row_ids, function(row_id) {
    .rc_step2_batch_target_spec(shared, payload, row_id)
  })
  names(specs) <- row_ids
  engine <- .rc_step2_batch_new_engine(shared, payload$solver)

  primary <- .rc_step2_batch_route_solve(
    specs = specs,
    engine = engine,
    penalty_matrix = payload$penalty,
    evidence = evidence,
    units = units
  )
  engine <- primary$engine
  control <- NULL
  if (paired_control) {
    control <- .rc_step2_batch_route_solve(
      specs = specs,
      engine = engine,
      penalty_matrix = payload$control_penalty,
      evidence = control_evidence,
      units = units,
      reuse_mask = reuse_mask,
      reuse_results = primary$results
    )
    engine <- control$engine
  }

  list(
    engine = engine,
    primary = primary,
    control = control,
    reuse_mask = reuse_mask,
    units = units,
    evidence = evidence,
    control_evidence = control_evidence
  )
}

.rc_step2_batch_primary_diagnostics <- function(
    row_id, entry, model, primary, routes, full_gem) {
  units <- routes$units
  if (isTRUE(full_gem)) {
    return(data.frame(
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
      step2_solver_reused_across_targets = rep(TRUE, length(units)),
      step2_traversal = rep("unit_then_directional_target", length(units)),
      target_expression_available = primary$target_available,
      objective_evidence_fraction = routes$evidence$fraction,
      unavailable_objective_terms = routes$evidence$unavailable,
      parallel_task = rep(
        "directional_reaction_x_all_metacells", length(units)
      ),
      stringsAsFactors = FALSE
    ))
  }

  target_status <- if (!is.null(model$target_status)) {
    rep(as.character(model$target_status), length(units))
  } else {
    ifelse(primary$feasible, "ok", "structurally_infeasible")
  }
  data.frame(
    row_id = rep(row_id, length(units)),
    unit_id = units,
    module_id = rep("CELLTYPE_MEDIUM_UNION_GEM", length(units)),
    cell_type = rep(entry$cell_type, length(units)),
    reaction_id = rep(entry$reaction_id, length(units)),
    target_direction = rep(entry$target_direction, length(units)),
    medium_scenario = rep(entry$medium_scenario, length(units)),
    condition = rep("all", length(units)),
    strict_feasible = unname(primary$feasible),
    solver_status = primary$solver_status,
    solver_backend = primary$solver_backend,
    step1_status = primary$step1_status,
    step2_status = primary$step2_status,
    target_status = target_status,
    objective_value = unname(primary$penalty),
    vmax = unname(primary$vmax),
    vmax_reused_from_celltype_cache = rep(TRUE, length(units)),
    step2_model_reused_across_metacells = rep(TRUE, length(units)),
    step2_solver_reused_across_targets = rep(TRUE, length(units)),
    step2_traversal = rep("unit_then_directional_target", length(units)),
    target_expression_available = primary$target_available,
    objective_evidence_fraction = routes$evidence$fraction,
    unavailable_objective_terms = routes$evidence$unavailable,
    parallel_task = rep(
      "directional_reaction_x_matching_metacells", length(units)
    ),
    stringsAsFactors = FALSE
  )
}

.rc_step2_batch_control_diagnostics <- function(
    row_id, entry, control, routes, full_gem) {
  if (is.null(control)) return(data.frame())
  units <- routes$units
  common <- list(
    row_id = rep(row_id, length(units)),
    unit_id = units,
    reaction_id = rep(entry$reaction_id, length(units)),
    target_direction = rep(entry$target_direction, length(units)),
    medium_scenario = rep(entry$medium_scenario, length(units)),
    objective_value = unname(control$penalty),
    strict_feasible = unname(control$feasible),
    solver_status = control$solver_status,
    solver_backend = control$solver_backend,
    objective_identical_to_primary = unname(routes$reuse_mask)
  )
  if (!isTRUE(full_gem)) {
    common <- append(
      common,
      list(cell_type = rep(entry$cell_type, length(units))),
      after = 2L
    )
  }
  as.data.frame(common, stringsAsFactors = FALSE)
}

.rc_step2_batch_checkpoint_worker <- function(task, full_gem = FALSE) {
  state <- .rc_step2_batch_validate_worker_payload(task, full_gem = full_gem)
  payload <- state$payload
  row_ids <- state$row_ids
  model <- state$model
  paired_control <- state$paired_control
  routes <- NULL
  step2_engine <- NULL
  on.exit({
    if (is.list(routes)) step2_engine <- routes$engine
    .rc_compass_step2_release_engine(step2_engine)
    rm(model, payload)
    invisible(gc(verbose = FALSE, full = TRUE))
  }, add = TRUE)

  routes <- .rc_step2_batch_run_routes(payload, row_ids, paired_control)
  step2_engine <- routes$engine
  units <- routes$units
  checkpoint_files <- character(length(row_ids))

  for (j in seq_along(row_ids)) {
    row_id <- row_ids[[j]]
    entry <- payload$entries[[row_id]]
    primary <- routes$primary$results[[row_id]]
    primary_metrics <- routes$primary$metrics[[row_id]]
    control <- if (paired_control) routes$control$results[[row_id]] else NULL
    control_metrics <- if (paired_control) {
      routes$control$metrics[[row_id]]
    } else {
      list(
        engine = "not_run", n_solves = 0L,
        n_objective_updates = 0L, n_fallback = 0L
      )
    }

    diagnostics <- .rc_step2_batch_primary_diagnostics(
      row_id, entry, model, primary, routes, full_gem
    )
    control_diagnostics <- .rc_step2_batch_control_diagnostics(
      row_id, entry, control, routes, full_gem
    )

    token <- substr(.rc_microcompass_object_checksum(list(
      row_id = row_id,
      file_checksum = model$file_checksum %||% NA_character_,
      units = units,
      omega = payload$omega,
      solver = payload$solver,
      flux_threshold = payload$flux_threshold
    )), 1L, 24L)
    checkpoint <- file.path(
      task$checkpoint_dir, paste0("step2__", token, ".rds")
    )

    primary_metrics$shared_model_batch_engine <- TRUE
    primary_metrics$batch_objective_change_events <-
      routes$primary$batch_metrics$n_objective_change_events
    primary_metrics$batch_target_switches <-
      routes$primary$batch_metrics$n_target_switches
    if (paired_control) {
      control_metrics$shared_model_batch_engine <- TRUE
      control_metrics$batch_objective_change_events <-
        routes$control$batch_metrics$n_objective_change_events
      control_metrics$batch_target_switches <-
        routes$control$batch_metrics$n_target_switches
    }

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
        reused_from_primary = isTRUE(all(routes$reuse_mask)),
        reused_from_primary_by_unit = routes$reuse_mask,
        shared_target_engine = TRUE,
        shared_model_batch_engine = TRUE
      ) else NULL
    ), checkpoint)
    checkpoint_files[[j]] <- checkpoint

    rm(
      diagnostics, control_diagnostics, primary, control,
      primary_metrics, control_metrics
    )
  }

  checkpoint_files
}

.rc_step2_reaction_batch_worker_model_batch <- function(task) {
  .rc_step2_batch_checkpoint_worker(task, full_gem = FALSE)
}

.rc_full_gem_step2_reaction_batch_worker_model_batch <- function(task) {
  .rc_step2_batch_checkpoint_worker(task, full_gem = TRUE)
}

# Production bindings. The established private worker symbols remain stable for
# the Layer 2 orchestration code, while each new implementation has one unique
# top-level function definition for API-surface auditing.
.rc_step2_reaction_batch_worker <- .rc_step2_reaction_batch_worker_model_batch
.rc_full_gem_step2_reaction_batch_worker <-
  .rc_full_gem_step2_reaction_batch_worker_model_batch
