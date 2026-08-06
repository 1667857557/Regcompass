# Original MATLAB CORDA2 three-step reconstruction state machine.

.rc_corda2_empty_dependency_matrix <- function(rows, columns) {
  matrix(
    0L,
    nrow = length(rows),
    ncol = length(columns),
    dimnames = list(rows, columns)
  )
}

.rc_corda2_reaction_numeric_confidence <- function(
    split, directional_class, included_variables = character()) {
  level <- c(OT = 0L, NC = 1L, MC = 2L, HC = 3L)
  directional <- unname(level[directional_class])
  names(directional) <- names(directional_class)
  if (length(included_variables)) {
    directional[] <- 0L
    directional[included_variables] <- 3L
  }
  reactions <- split$reaction_order
  answer <- stats::setNames(integer(length(reactions)), reactions)
  for (reaction in reactions) {
    variables <- split$direction_table$variable_id[
      split$direction_table$reaction_id == reaction
    ]
    answer[[reaction]] <- max(directional[variables])
  }
  answer
}

.rc_corda_build_three_stage_core <- function(
    split, classes, options, solver, time_limit) {
  directional_class <- .rc_corda2_directional_class(split, classes)
  initial_directional_class <- directional_class
  hc <- names(directional_class)[directional_class == "HC"]
  mc <- names(directional_class)[directional_class == "MC"]
  nc <- names(directional_class)[directional_class == "NC"]
  ot <- names(directional_class)[directional_class == "OT"]
  inclusion_stage_direction <- stats::setNames(
    rep(NA_character_, length(directional_class)), names(directional_class)
  )
  inclusion_stage_direction[hc] <- "initial_high_confidence"

  engine <- .rc_corda_new_lp_engine(split, solver, time_limit)
  on.exit({
    engine <- .rc_corda_release_lp_engine(engine)
  }, add = TRUE)
  task_tables <- list()

  # Step 1: support every high-confidence direction with MC and NC reactions.
  stage1_hc <- hc
  stage1_mc <- mc
  stage1_nc <- nc
  HCtoMC <- .rc_corda2_empty_dependency_matrix(stage1_hc, stage1_mc)
  HCtoNC <- .rc_corda2_empty_dependency_matrix(stage1_hc, stage1_nc)
  hc_present <- mc_present <- nc_present <- character()
  blocked_hc <- logical(length(stage1_hc))
  stage1_results <- vector("list", length(stage1_hc))
  for (i in seq_along(stage1_hc)) {
    target <- stage1_hc[[i]]
    assessed <- .rc_corda2_dependency_assessment(
      engine, split, target, directional_class, options,
      stage = "corda2_step1_HC_dependencies",
      penalized_class = "stage1"
    )
    engine <- assessed$engine
    stage1_results[[i]] <- assessed$result
    if (!isTRUE(assessed$success)) {
      blocked_hc[[i]] <- TRUE
      next
    }
    active <- assessed$active
    hc_present <- union(hc_present, active[directional_class[active] == "HC"])
    mc_used <- assessed$associated[
      directional_class[assessed$associated] == "MC"
    ]
    nc_used <- assessed$associated[
      directional_class[assessed$associated] == "NC"
    ]
    mc_present <- union(mc_present, mc_used)
    nc_present <- union(nc_present, nc_used)
    if (length(mc_used)) HCtoMC[i, mc_used] <- 1L
    if (length(nc_used)) HCtoNC[i, nc_used] <- 1L
  }
  blocked_hc[stage1_hc %in% hc_present] <- FALSE
  retained_stage1_hc <- stage1_hc[!blocked_hc]
  HCtoMC <- HCtoMC[!blocked_hc, , drop = FALSE]
  HCtoNC <- HCtoNC[!blocked_hc, , drop = FALSE]
  promoted_step1_mc <- stage1_mc[stage1_mc %in% mc_present]
  promoted_step1_nc <- stage1_nc[stage1_nc %in% nc_present]
  promoted_step1 <- c(promoted_step1_mc, promoted_step1_nc)
  inclusion_stage_direction[promoted_step1] <-
    "corda2_step1_associated_with_HC"
  directional_class[promoted_step1] <- "HC"
  hc <- c(retained_stage1_hc, promoted_step1)
  mc <- stage1_mc[!stage1_mc %in% mc_present]
  nc <- stage1_nc[!stage1_nc %in% nc_present]
  task_tables$step1 <- .rc_corda2_results_table(stage1_results)

  # Step 2.1: determine NC dependencies of every remaining MC direction.
  stage2_mc_input <- mc
  stage2_nc_input <- nc
  MCxNC <- .rc_corda2_empty_dependency_matrix(stage2_mc_input, stage2_nc_input)
  blocked_mc_step21 <- logical(length(stage2_mc_input))
  stage21_results <- vector("list", length(stage2_mc_input))
  for (i in seq_along(stage2_mc_input)) {
    target <- stage2_mc_input[[i]]
    assessed <- .rc_corda2_dependency_assessment(
      engine, split, target, directional_class, options,
      stage = "corda2_step2_1_MC_NC_dependencies",
      penalized_class = "NC"
    )
    engine <- assessed$engine
    stage21_results[[i]] <- assessed$result
    if (!isTRUE(assessed$success)) {
      blocked_mc_step21[[i]] <- TRUE
      next
    }
    nc_used <- assessed$associated[
      directional_class[assessed$associated] == "NC"
    ]
    if (length(nc_used)) MCxNC[i, nc_used] <- 1L
  }
  MCxNC <- MCxNC[!blocked_mc_step21, , drop = FALSE]
  mc <- stage2_mc_input[!blocked_mc_step21]
  MCtoNC <- MCxNC
  task_tables$step2_1 <- .rc_corda2_results_table(stage21_results)

  # Step 2.2: promote frequently required NC directions, block the rest, and
  # retain only MC directions that remain feasible.
  nc_count <- if (ncol(MCxNC)) colSums(MCxNC) else numeric()
  promoted_nc <- names(nc_count)[nc_count >= options$MCxNCthresh]
  if (length(promoted_nc)) {
    promoted_rows <- matrix(
      0L, nrow = length(promoted_nc), ncol = ncol(MCxNC),
      dimnames = list(promoted_nc, colnames(MCxNC))
    )
    MCxNC <- rbind(MCxNC, promoted_rows)
  }
  directional_class[promoted_nc] <- "MC"
  mc <- c(mc, promoted_nc)
  if (length(promoted_nc) && ncol(MCxNC)) {
    MCxNC <- MCxNC[, setdiff(colnames(MCxNC), promoted_nc), drop = FALSE]
  }
  nc <- setdiff(nc, promoted_nc)
  split_step22 <- split
  if (length(nc)) {
    split_step22$lb[nc] <- 0
    split_step22$ub[nc] <- 0
  }
  stage22_results <- vector("list", length(mc))
  blocked_mc_step22 <- logical(length(mc))
  rescue <- vector("list", length(mc))
  for (i in seq_along(mc)) {
    target <- mc[[i]]
    maximum <- .rc_corda2_maximize_target(
      engine, split_step22, target,
      lower = split_step22$lb,
      upper = split_step22$ub
    )
    engine <- maximum$engine
    success <- identical(maximum$answer$status, "optimal") &&
      is.finite(maximum$vmax) && maximum$vmax >= options$flux_threshold
    blocked_mc_step22[[i]] <- !success
    dependencies <- if (target %in% rownames(MCxNC)) {
      colnames(MCxNC)[MCxNC[target, , drop = TRUE] > 0]
    } else {
      character()
    }
    rescue[[i]] <- dependencies
    stage22_results[[i]] <- .rc_corda2_target_result(
      split, target, "corda2_step2_2_MC_feasibility",
      if (success) "optimal" else "target_blocked",
      associated = dependencies,
      target_flux = maximum$vmax,
      vmax = maximum$vmax,
      objective = maximum$answer$objective,
      backend = maximum$answer$backend,
      solver_message = maximum$answer$solver_message %||% "",
      n_solves = 1L,
      opposite = maximum$opposite
    )
  }
  deleted_mc <- mc[blocked_mc_step22]
  rescue_table <- data.frame(
    reaction = deleted_mc,
    dependent_on = vapply(
      rescue[blocked_mc_step22], paste, character(1), collapse = ","
    ),
    stringsAsFactors = FALSE
  )
  feasible_mc <- mc[!blocked_mc_step22]
  inclusion_stage_direction[promoted_nc] <-
    "corda2_step2_2_NC_occurrence_threshold"
  inclusion_stage_direction[feasible_mc[is.na(
    inclusion_stage_direction[feasible_mc]
  )]] <- "corda2_step2_2_MC_feasible"
  directional_class[feasible_mc] <- "HC"
  hc <- c(hc, feasible_mc)
  task_tables$step2_2 <- .rc_corda2_results_table(stage22_results)

  # Step 3: block all remaining MC/NC directions and add only OT reactions
  # required for retained HC flux.
  split_step3 <- split_step22
  allowed_step3 <- union(hc, ot)
  blocked_step3 <- setdiff(colnames(split_step3$S), allowed_step3)
  if (length(blocked_step3)) {
    split_step3$lb[blocked_step3] <- 0
    split_step3$ub[blocked_step3] <- 0
  }
  stage3_results <- vector("list", length(hc))
  ot_present <- character()
  for (i in seq_along(hc)) {
    target <- hc[[i]]
    assessed <- .rc_corda2_dependency_assessment(
      engine, split_step3, target, directional_class, options,
      stage = "corda2_step3_HC_OT_dependencies",
      penalized_class = "OT",
      lower = split_step3$lb,
      upper = split_step3$ub
    )
    engine <- assessed$engine
    stage3_results[[i]] <- assessed$result
    if (!isTRUE(assessed$success)) next
    used <- assessed$associated[
      directional_class[assessed$associated] == "OT"
    ]
    ot_present <- union(ot_present, used)
  }
  inclusion_stage_direction[ot_present] <-
    "corda2_step3_associated_OT"
  included_variables <- unique(c(hc, ot_present))
  task_tables$step3 <- .rc_corda2_results_table(stage3_results)

  initial_reaction_confidence <- .rc_corda2_reaction_numeric_confidence(
    split, initial_directional_class
  )
  final_reaction_confidence <- .rc_corda2_reaction_numeric_confidence(
    split, initial_directional_class, included_variables
  )
  included_reactions <- unique(as.character(
    split$variable_to_reaction[included_variables]
  ))
  final_reaction_status <- stats::setNames(
    ifelse(names(final_reaction_confidence) %in% included_reactions,
           "included", "excluded"),
    names(final_reaction_confidence)
  )
  inclusion_stage <- stats::setNames(
    rep(NA_character_, length(split$reaction_order)), split$reaction_order
  )
  for (reaction in names(inclusion_stage)) {
    variables <- split$direction_table$variable_id[
      split$direction_table$reaction_id == reaction &
        split$direction_table$variable_id %in% included_variables
    ]
    stages <- inclusion_stage_direction[variables]
    stages <- stages[!is.na(stages)]
    if (length(stages)) inclusion_stage[[reaction]] <- stages[[1L]]
  }

  list(
    included = included_reactions,
    included_directional_variables = included_variables,
    initial_reaction_confidence = initial_reaction_confidence,
    final_reaction_confidence = final_reaction_confidence,
    final_reaction_status = final_reaction_status,
    final_confidence = final_reaction_confidence,
    initial_directional_class = initial_directional_class,
    inclusion_stage = inclusion_stage,
    inclusion_stage_direction = inclusion_stage_direction,
    stage1_associated = unique(as.character(
      split$variable_to_reaction[promoted_step1]
    )),
    stage2_promoted_nc = unique(as.character(
      split$variable_to_reaction[promoted_nc]
    )),
    stage2_promoted_mc = unique(as.character(
      split$variable_to_reaction[feasible_mc]
    )),
    stage3_associated_ot = unique(as.character(
      split$variable_to_reaction[ot_present]
    )),
    blocked_after_stage2 = unique(as.character(
      split$variable_to_reaction[nc]
    )),
    blocked_before_stage3 = unique(as.character(
      split$variable_to_reaction[blocked_step3]
    )),
    blocked_high_confidence_directions = stage1_hc[blocked_hc],
    blocked_medium_directions_step2_1 = stage2_mc_input[blocked_mc_step21],
    blocked_medium_directions_step2_2 = deleted_mc,
    HCtoMC = HCtoMC,
    HCtoNC = HCtoNC,
    MCtoNC = MCtoNC,
    rescue = rescue_table,
    task_diagnostics = .rc_bind_frames_fill(task_tables),
    solver_performance = .rc_corda_execution_metrics(engine),
    algorithm = "schultzdre_MATLAB_CORDA2_original_semantics",
    reference_repository = "schultzdre/Constraint-Based-Modeling",
    reference_file = "CORDA2.m",
    stage_update_policy = "original_matlab_directional_order",
    source_semantics = c(
      "split only actively reversible reactions once",
      "close the opposite direction for every tested reaction",
      "fix the target at the original val-or-percentage constraint",
      "step 1 costs MC=sqrt(om), NC=om and all other reactions=1e-3",
      "increase each newly used high-cost reaction once by 1+ci",
      "promote NC directions required by at least MCxNCthresh MC directions",
      "block remaining NC before MC feasibility and add OT only in step 3"
    )
  )
}

# Progress-aware entry point; the algorithm remains in the core above.
.rc_corda_build_three_stage <- function(...) {
  answer <- do.call(
    .rc_corda_build_three_stage_core,
    list(...)
  )
  task <- .rc_layer2_progress_state$current_task
  if (!is.null(task) && identical(task$route, "corda2")) {
    required <- list(
      c("corda2_step1", "corda2_step1_HC_dependencies", 4L),
      c("corda2_step2_1", "corda2_step2_1_MC_NC_dependencies", 5L),
      c("corda2_step2_2", "corda2_step2_2_MC_feasibility", 6L),
      c("corda2_step3", "corda2_step3_HC_OT_dependencies", 7L)
    )
    for (item in required) {
      if (!exists(item[[1L]], envir = .rc_layer2_progress_state$algorithm_flags,
                  inherits = FALSE)) {
        .rc_layer2_algorithm_once(
          item[[1L]], item[[2L]], as.integer(item[[3L]]),
          "skipped because this confidence class had no candidate directions"
        )
      }
    }
  }
  answer
}
