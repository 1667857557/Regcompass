# Bounded-memory production wrapper for model-batch Layer 2 COMPASS Step 2.
#
# The exact LP formulation and target-bound logic remain in the canonical
# model-batch helpers. This wrapper keeps one solver for the whole worker, but
# limits the number of target-by-unit result vectors retained at once. It also
# reports cross-target persistent-solver reuse only when that reuse actually
# occurred. CORDA2, directional Vmax, penalties, omega, media and targets are
# unchanged.

.rc_step2_stream_target_chunk_size <- function(n_units, n_targets) {
  n_units <- max(1L, as.integer(n_units[[1L]]))
  n_targets <- max(0L, as.integer(n_targets[[1L]]))
  if (!n_targets) return(0L)

  # Each route result stores several vectors of length n_units per target.
  # Keep the additional target-by-unit state bounded independently of the
  # number of targets assigned to a worker. A single target is the irreducible
  # lower bound when n_units itself exceeds this engineering budget.
  max_target_unit_pairs <- 100000L
  max_targets <- 16L
  by_pairs <- max(1L, floor(max_target_unit_pairs / n_units))
  as.integer(min(n_targets, max_targets, by_pairs))
}

.rc_step2_stream_cross_target_reuse <- function(
    engine_type, runnable_target_count, row_runnable, n_solves) {
  identical(as.character(engine_type %||% "not_run"),
            "highs_persistent_cpp") &&
    as.integer(runnable_target_count %||% 0L) > 1L &&
    isTRUE(row_runnable) &&
    as.integer(n_solves %||% 0L) > 0L
}

.rc_step2_stream_primary_diagnostics <- function(
    row_id, entry, model, primary, routes, full_gem, solver_reused) {
  out <- .rc_step2_batch_primary_diagnostics(
    row_id, entry, model, primary, routes, full_gem
  )
  out$step2_solver_reused_across_targets <-
    rep(isTRUE(solver_reused), nrow(out))
  out$step2_traversal <- rep(
    "bounded_target_microbatch_unit_then_directional_target", nrow(out)
  )
  out
}

.rc_step2_stream_checkpoint_worker <- function(task, full_gem = FALSE) {
  state <- .rc_step2_batch_validate_worker_payload(task, full_gem = full_gem)
  payload <- state$payload
  row_ids <- state$row_ids
  model <- state$model
  paired_control <- state$paired_control
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

  # Compile the target-invariant directional model exactly once for the full
  # worker batch. Target specifications are O(number of targets) small metadata;
  # the large O(targets x units) result vectors are never retained globally.
  shared <- .rc_step2_batch_compile_shared(payload, row_ids)
  specs <- lapply(row_ids, function(row_id) {
    .rc_step2_batch_target_spec(shared, payload, row_id)
  })
  names(specs) <- row_ids
  runnable_target_count <- sum(vapply(
    specs, function(spec) isTRUE(spec$runnable), logical(1)
  ))
  step2_engine <- .rc_step2_batch_new_engine(shared, payload$solver)
  engine_type <- if (is.list(step2_engine)) {
    as.character(step2_engine$type %||% "one_shot")
  } else {
    "not_run"
  }

  on.exit({
    .rc_compass_step2_release_engine(step2_engine)
    rm(model, payload, shared, specs)
    invisible(gc(verbose = FALSE, full = TRUE))
  }, add = TRUE)

  target_chunk_size <- .rc_step2_stream_target_chunk_size(
    length(units), length(row_ids)
  )
  chunk_id <- ceiling(seq_along(row_ids) / target_chunk_size)
  chunks <- unname(split(row_ids, chunk_id))
  checkpoint_files <- stats::setNames(character(length(row_ids)), row_ids)
  worker_primary_objective_events <- 0L
  worker_primary_target_switches <- 0L
  worker_control_objective_events <- 0L
  worker_control_target_switches <- 0L

  for (chunk_index in seq_along(chunks)) {
    chunk_rows <- as.character(chunks[[chunk_index]])
    chunk_specs <- specs[chunk_rows]

    primary <- .rc_step2_batch_route_solve(
      specs = chunk_specs,
      engine = step2_engine,
      penalty_matrix = payload$penalty,
      evidence = evidence,
      units = units
    )
    step2_engine <- primary$engine
    worker_primary_objective_events <-
      worker_primary_objective_events +
      as.integer(primary$batch_metrics$n_objective_change_events %||% 0L)
    worker_primary_target_switches <-
      worker_primary_target_switches +
      as.integer(primary$batch_metrics$n_target_switches %||% 0L)

    control <- NULL
    if (paired_control) {
      control <- .rc_step2_batch_route_solve(
        specs = chunk_specs,
        engine = step2_engine,
        penalty_matrix = payload$control_penalty,
        evidence = control_evidence,
        units = units,
        reuse_mask = reuse_mask,
        reuse_results = primary$results
      )
      step2_engine <- control$engine
      worker_control_objective_events <-
        worker_control_objective_events +
        as.integer(control$batch_metrics$n_objective_change_events %||% 0L)
      worker_control_target_switches <-
        worker_control_target_switches +
        as.integer(control$batch_metrics$n_target_switches %||% 0L)
    }

    routes <- list(
      engine = step2_engine,
      primary = primary,
      control = control,
      reuse_mask = reuse_mask,
      units = units,
      evidence = evidence,
      control_evidence = control_evidence
    )

    # Checkpoint every row as soon as its bounded micro-batch has finished.
    # Only this chunk's primary/control result vectors remain resident.
    for (row_id in chunk_rows) {
      entry <- payload$entries[[row_id]]
      spec <- chunk_specs[[row_id]]
      primary_result <- primary$results[[row_id]]
      primary_metrics <- primary$metrics[[row_id]]
      control_result <- if (paired_control) control$results[[row_id]] else NULL
      control_metrics <- if (paired_control) {
        control$metrics[[row_id]]
      } else {
        list(
          engine = "not_run", n_solves = 0L,
          n_objective_updates = 0L, n_fallback = 0L,
          n_bound_updates = 0L, n_target_switches = 0L
        )
      }

      primary_shared <- .rc_step2_stream_cross_target_reuse(
        primary_metrics$engine,
        runnable_target_count,
        spec$runnable,
        primary_metrics$n_solves
      )
      control_shared <- if (paired_control) {
        .rc_step2_stream_cross_target_reuse(
          control_metrics$engine,
          runnable_target_count,
          spec$runnable,
          control_metrics$n_solves
        )
      } else {
        FALSE
      }

      diagnostics <- .rc_step2_stream_primary_diagnostics(
        row_id, entry, model, primary_result, routes, full_gem,
        solver_reused = primary_shared
      )
      control_diagnostics <- .rc_step2_batch_control_diagnostics(
        row_id, entry, control_result, routes, full_gem
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

      primary_metrics$shared_model_batch_engine <- isTRUE(primary_shared)
      primary_metrics$batch_objective_change_events <-
        as.integer(primary$batch_metrics$n_objective_change_events %||% 0L)
      primary_metrics$batch_target_switches <-
        as.integer(primary$batch_metrics$n_target_switches %||% 0L)
      primary_metrics$worker_objective_change_events <-
        as.integer(worker_primary_objective_events)
      primary_metrics$worker_target_switches <-
        as.integer(worker_primary_target_switches)
      primary_metrics$stream_target_chunk_size <-
        as.integer(target_chunk_size)
      primary_metrics$stream_target_chunk_count <-
        as.integer(length(chunks))
      primary_metrics$stream_chunk_target_unit_pairs <-
        as.double(length(chunk_rows)) * as.double(length(units))

      if (paired_control) {
        control_metrics$shared_model_batch_engine <- isTRUE(control_shared)
        control_metrics$batch_objective_change_events <-
          as.integer(control$batch_metrics$n_objective_change_events %||% 0L)
        control_metrics$batch_target_switches <-
          as.integer(control$batch_metrics$n_target_switches %||% 0L)
        control_metrics$worker_objective_change_events <-
          as.integer(worker_control_objective_events)
        control_metrics$worker_target_switches <-
          as.integer(worker_control_target_switches)
        control_metrics$stream_target_chunk_size <-
          as.integer(target_chunk_size)
        control_metrics$stream_target_chunk_count <-
          as.integer(length(chunks))
        control_metrics$stream_chunk_target_unit_pairs <-
          as.double(length(chunk_rows)) * as.double(length(units))
      }

      paired_shared_target_engine <- paired_control &&
        identical(engine_type, "highs_persistent_cpp") &&
        isTRUE(spec$runnable) &&
        as.integer(primary_metrics$n_solves %||% 0L) > 0L &&
        as.integer(control_metrics$n_solves %||% 0L) > 0L

      .rc_atomic_save_rds(list(
        row_id = row_id,
        units = units,
        penalty = primary_result$penalty,
        vmax = primary_result$vmax,
        feasible = primary_result$feasible,
        evaluated = primary_result$evaluated,
        diagnostics = diagnostics,
        engine_metrics = primary_metrics,
        control = if (paired_control) list(
          penalty = control_result$penalty,
          vmax = control_result$vmax,
          feasible = control_result$feasible,
          evaluated = control_result$evaluated,
          diagnostics = control_diagnostics,
          engine_metrics = control_metrics,
          reused_from_primary = isTRUE(all(reuse_mask)),
          reused_from_primary_by_unit = reuse_mask,
          shared_target_engine = isTRUE(paired_shared_target_engine),
          shared_model_batch_engine = isTRUE(control_shared)
        ) else NULL
      ), checkpoint)
      checkpoint_files[[row_id]] <- checkpoint

      # Release completed row vectors immediately within the micro-batch.
      primary$results[[row_id]] <- NULL
      primary$metrics[[row_id]] <- NULL
      if (paired_control) {
        control$results[[row_id]] <- NULL
        control$metrics[[row_id]] <- NULL
      }
      rm(
        diagnostics, control_diagnostics, primary_result, control_result,
        primary_metrics, control_metrics
      )
    }

    rm(routes, primary, control, chunk_specs)
    invisible(gc(verbose = FALSE, full = FALSE))
  }

  unname(checkpoint_files[row_ids])
}

.rc_step2_reaction_batch_worker_streaming <- function(task) {
  .rc_step2_stream_checkpoint_worker(task, full_gem = FALSE)
}

.rc_full_gem_step2_reaction_batch_worker_streaming <- function(task) {
  .rc_step2_stream_checkpoint_worker(task, full_gem = TRUE)
}

# Final production bindings are collated after layer2_step2_model_batch.R.
.rc_step2_reaction_batch_worker <- .rc_step2_reaction_batch_worker_streaming
.rc_full_gem_step2_reaction_batch_worker <-
  .rc_full_gem_step2_reaction_batch_worker_streaming
