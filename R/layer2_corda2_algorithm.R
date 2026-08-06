# Original MATLAB CORDA2 dependency assessment.

.rc_corda2_directional_class <- function(split, classes) {
  reaction_class <- stats::setNames(
    rep("OT", length(classes$confidence)), names(classes$confidence)
  )
  reaction_class[classes$hc] <- "HC"
  reaction_class[classes$mc] <- "MC"
  reaction_class[classes$nc] <- "NC"
  reaction_class[classes$ot] <- "OT"
  value <- unname(reaction_class[split$direction_table$reaction_id])
  if (anyNA(value) || any(!value %in% c("HC", "MC", "NC", "OT"))) {
    stop("Every CORDA2 directional variable requires HC, MC, NC or OT class.",
         call. = FALSE)
  }
  stats::setNames(value, split$direction_table$variable_id)
}

.rc_corda2_stage_cost <- function(
    split, directional_class, options,
    penalized_class, baseline = 0) {
  cost <- stats::setNames(rep(baseline, ncol(split$S)), colnames(split$S))
  if (identical(penalized_class, "stage1")) {
    cost[directional_class == "MC"] <- sqrt(options$om)
    cost[directional_class == "NC"] <- options$om
  } else if (identical(penalized_class, "NC")) {
    cost[directional_class == "NC"] <- options$om
  } else if (identical(penalized_class, "OT")) {
    cost[directional_class == "OT"] <- options$om
  } else {
    stop("Unknown CORDA2 cost stage.", call. = FALSE)
  }
  cost
}

.rc_corda2_target_result <- function(
    split, target, stage, status, associated = character(),
    active = character(), target_flux = NA_real_, vmax = NA_real_,
    objective = NA_real_, backend = "", solver_message = "",
    n_solves = 0L, opposite = character()) {
  row <- split$direction_table[
    split$direction_table$variable_id == target, , drop = FALSE
  ]
  list(
    target = target,
    reaction_id = as.character(row$reaction_id[[1L]]),
    direction = as.character(row$direction[[1L]]),
    stage = stage,
    kind = "dependency",
    status = status,
    associated = as.character(associated),
    active = as.character(active),
    target_flux = as.numeric(target_flux),
    vmax = as.numeric(vmax),
    objective = as.numeric(objective),
    backend = as.character(backend),
    solver_message = as.character(solver_message),
    opposite_direction_blocked = as.character(opposite),
    redundancies = max(0L, as.integer(n_solves) - 1L),
    n_solves = as.integer(n_solves)
  )
}

.rc_corda2_scan_flux <- function(
    flux, class_code, track_code, threshold) {
  native <- get0(
    ".rc_corda2_scan_flux_cpp", mode = "function", inherits = TRUE
  )
  if (is.function(native)) {
    return(native(flux, class_code, track_code, threshold))
  }
  active <- which(is.finite(flux) & flux > threshold)
  used <- active[class_code[active] %in% track_code]
  list(active = as.integer(active), used = as.integer(used))
}

.rc_corda2_dependency_assessment_core <- function(
    engine, split, target, directional_class, options,
    stage, penalized_class,
    lower = split$lb, upper = split$ub,
    constrain_target = TRUE) {
  constrained <- if (isTRUE(constrain_target)) {
    .rc_corda2_constrain_target(
      engine, split, target, options, lower = lower, upper = upper
    )
  } else {
    .rc_corda2_maximize_target(
      engine, split, target, lower = lower, upper = upper
    )
  }
  engine <- constrained$engine
  if (!identical(constrained$answer$status, "optimal") ||
      !is.finite(constrained$vmax) ||
      constrained$vmax < options$flux_threshold ||
      (isTRUE(constrain_target) && !is.finite(constrained$required_flux))) {
    return(list(
      engine = engine,
      result = .rc_corda2_target_result(
        split, target, stage, "target_blocked",
        target_flux = 0,
        vmax = constrained$vmax,
        backend = constrained$answer$backend,
        solver_message = constrained$answer$solver_message %||% "",
        n_solves = 1L,
        opposite = constrained$opposite
      ),
      associated = character(),
      active = character(),
      success = FALSE
    ))
  }

  baseline <- if (identical(penalized_class, "stage1")) {
    options$baseline_cost
  } else {
    0
  }
  penalty <- .rc_corda2_stage_cost(
    split, directional_class, options,
    penalized_class = penalized_class,
    baseline = baseline
  )
  track_class <- switch(
    penalized_class,
    stage1 = c("MC", "NC"),
    NC = "NC",
    OT = "OT"
  )
  variable_ids <- colnames(split$S)
  target_index <- match(target, variable_ids)
  class_levels <- c("HC", "MC", "NC", "OT")
  class_code <- match(unname(directional_class[variable_ids]), class_levels)
  track_code <- match(track_class, class_levels)
  if (is.na(target_index) || anyNA(class_code) || anyNA(track_code)) {
    stop("CORDA2 native dependency indices are incomplete.", call. = FALSE)
  }

  associated <- character()
  active_all <- character()
  associated_seen <- rep(FALSE, length(variable_ids))
  active_seen <- rep(FALSE, length(variable_ids))
  n_solves <- 1L
  final_answer <- constrained$answer
  final_flux <- NA_real_

  repeat {
    solved <- .rc_corda_engine_solve(
      engine,
      objective = as.numeric(penalty),
      lower = constrained$lower,
      upper = constrained$upper
    )
    engine <- solved$engine
    answer <- solved$answer
    n_solves <- n_solves + 1L
    final_answer <- answer
    if (!identical(answer$status, "optimal") ||
        length(answer$solution) != ncol(split$S)) {
      return(list(
        engine = engine,
        result = .rc_corda2_target_result(
          split, target, stage, answer$status,
          associated = associated,
          active = active_all,
          vmax = constrained$vmax,
          objective = answer$objective,
          backend = answer$backend,
          solver_message = answer$solver_message %||% "",
          n_solves = n_solves,
          opposite = constrained$opposite
        ),
        associated = associated,
        active = active_all,
        success = FALSE
      ))
    }
    flux <- as.numeric(answer$solution)
    final_flux <- flux[[target_index]]
    if (!is.finite(final_flux) || final_flux < options$flux_threshold) {
      return(list(
        engine = engine,
        result = .rc_corda2_target_result(
          split, target, stage, "target_below_flux_threshold",
          associated = associated,
          active = active_all,
          target_flux = final_flux,
          vmax = constrained$vmax,
          objective = answer$objective,
          backend = answer$backend,
          solver_message = answer$solver_message %||% "",
          n_solves = n_solves,
          opposite = constrained$opposite
        ),
        associated = associated,
        active = active_all,
        success = FALSE
      ))
    }

    scan <- .rc_corda2_scan_flux(
      flux = flux,
      class_code = class_code,
      track_code = track_code,
      threshold = options$flux_threshold
    )
    active_index <- as.integer(scan$active)
    used_index <- as.integer(scan$used)

    new_active <- active_index[!active_seen[active_index]]
    if (length(new_active)) {
      active_seen[new_active] <- TRUE
      active_all <- c(active_all, variable_ids[new_active])
    }

    newly_used <- used_index[!associated_seen[used_index]]
    if (length(newly_used)) {
      associated_seen[newly_used] <- TRUE
      associated <- c(associated, variable_ids[newly_used])
      penalty[newly_used] <- penalty[newly_used] * (1 + options$ci)
    } else {
      break
    }
  }

  list(
    engine = engine,
    result = .rc_corda2_target_result(
      split, target, stage, final_answer$status,
      associated = associated,
      active = active_all,
      target_flux = final_flux,
      vmax = constrained$vmax,
      objective = final_answer$objective,
      backend = final_answer$backend,
      solver_message = final_answer$solver_message %||% "",
      n_solves = n_solves,
      opposite = constrained$opposite
    ),
    associated = associated,
    active = active_all,
    success = TRUE
  )
}

.rc_corda2_results_table <- function(results, split = NULL) {
  if (!length(results)) return(data.frame())
  rows <- lapply(results, function(x) {
    data.frame(
      variable_id = as.character(x$target),
      reaction_id = as.character(x$reaction_id),
      direction = as.character(x$direction),
      stage = as.character(x$stage),
      replicate = 1L,
      kind = as.character(x$kind),
      status = as.character(x$status),
      target_flux = as.numeric(x$target_flux),
      vmax = as.numeric(x$vmax),
      objective = as.numeric(x$objective),
      backend = as.character(x$backend),
      solver_message = as.character(x$solver_message %||% ""),
      opposite_direction_blocked = paste(
        x$opposite_direction_blocked, collapse = ";"
      ),
      n_associated = length(x$associated),
      associated = paste(x$associated, collapse = ";"),
      corda2_n_solves = as.integer(x$n_solves %||% 0L),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Progress-aware entry point; the algorithm remains in the core above.
.rc_corda2_dependency_assessment <- function(...) {
  progress_state <- get0(
    ".rc_layer2_progress_state", mode = "environment", inherits = TRUE
  )
  args <- list(...)
  stage <- args$stage %||% if (length(args) >= 6L) args[[6L]] else ""
  task <- progress_state$current_task
  if (!is.null(task) && identical(task$route, "corda2")) {
    if (identical(stage, "corda2_step1_HC_dependencies")) {
      .rc_layer2_algorithm_once(
        "corda2_step1", "corda2_step1_HC_dependencies", 4L,
        "supporting high-confidence directions with MC/NC dependencies"
      )
    } else if (identical(stage, "corda2_step2_1_MC_NC_dependencies")) {
      .rc_layer2_algorithm_once(
        "corda2_step2_1", "corda2_step2_1_MC_NC_dependencies", 5L,
        "measuring NC dependencies of remaining MC directions"
      )
    } else if (identical(stage, "corda2_step3_HC_OT_dependencies")) {
      .rc_layer2_algorithm_once(
        "corda2_step3", "corda2_step3_HC_OT_dependencies", 7L,
        "adding only OT reactions required by retained HC flux"
      )
    }
  }
  previous <- progress_state$inside_dependency
  progress_state$inside_dependency <- TRUE
  on.exit({ progress_state$inside_dependency <- previous }, add = TRUE)
  do.call(
    .rc_corda2_dependency_assessment_core,
    args
  )
}
