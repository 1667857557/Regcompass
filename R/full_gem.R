.rc_validate_full_gem_model_params <- function(model_params = list()) {
  if (!is.list(model_params)) {
    stop("`model_params` must be a list.", call. = FALSE)
  }
  requested <- as.character(model_params$model_completion %||% "none")
  if (length(requested) != 1L || is.na(requested) ||
      !identical(requested, "none")) {
    stop(
      "Full-GEM mode automatically skips FASTCORE and CORDA2; omit ",
      "`model_completion` or set it to `none`.",
      call. = FALSE
    )
  }
  incompatible <- intersect(names(model_params), c(
    "fastcore_epsilon", "max_support_reactions", "strict",
    "corda2_args", "corda_medium_confidence_threshold",
    "corda_negative_confidence_threshold", "corda_regulatory_weight",
    "corda_include_evidence_outside_modules",
    "corda_max_medium_confidence_reactions"
  ))
  if (length(incompatible)) {
    stop(
      "Full-GEM mode does not accept FASTCORE or CORDA2 controls: ",
      paste(incompatible, collapse = ", "),
      ". It applies medium bounds to the complete GEM and evaluates ",
      "directional feasibility during COMPASS scoring.",
      call. = FALSE
    )
  }
  unknown <- setdiff(
    names(model_params),
    c("model_completion", "completion_time_limit", "cache_dir")
  )
  if (length(unknown)) {
    stop(
      "Unsupported full-GEM `model_params`: ",
      paste(unknown, collapse = ", "),
      ".", call. = FALSE
    )
  }
  if (!is.null(model_params$completion_time_limit)) {
    value <- as.numeric(model_params$completion_time_limit)
    if (length(value) != 1L || !is.finite(value) || value <= 0) {
      stop("`completion_time_limit` must be positive.", call. = FALSE)
    }
  }
  invisible(list(model_completion = "none"))
}

.rc_assert_medium_bounds_only <- function(
    reference_gem, constrained_gem, medium_table = NULL,
    context = "medium application") {
  reference <- rc_validate_gem(reference_gem)
  constrained <- rc_validate_gem(constrained_gem)

  same_structure <- identical(reference$reactions, constrained$reactions) &&
    identical(reference$metabolites, constrained$metabolites) &&
    identical(reference$S, constrained$S)
  if (!same_structure) {
    stop(
      context,
      " changed the reaction/metabolite set or stoichiometric matrix. ",
      "Medium handling must only modify reaction bounds.",
      call. = FALSE
    )
  }

  changed <- reference$reactions[
    reference$lb != constrained$lb | reference$ub != constrained$ub
  ]
  allowed <- if (!is.null(medium_table) &&
                 "exchange_reaction_id" %in% colnames(medium_table)) {
    unique(as.character(medium_table$exchange_reaction_id))
  } else {
    character()
  }
  allowed <- allowed[!is.na(allowed) & nzchar(allowed)]
  unexpected <- setdiff(changed, allowed)
  if (length(unexpected)) {
    stop(
      context,
      " modified bounds outside the supplied exchange reactions: ",
      paste(utils::head(unexpected, 10L), collapse = ", "),
      ".", call. = FALSE
    )
  }

  invisible(list(
    changed_reactions = changed,
    n_changed_reactions = length(changed),
    n_removed_reactions = 0L
  ))
}

#' Build a COMPASS-style full GEM for one medium scenario
#'
#' The complete validated GEM is retained. Optional medium constraints modify
#' exchange-reaction bounds only. Reactions that cannot carry flux under the
#' medium remain in the model and are classified by the directional COMPASS
#' maximum-flux step; the penalty-minimization step is skipped when infeasible.
.rc_build_full_gem_core <- function(gem, medium_table = NULL, condition = NULL) {
  gem <- rc_annotate_reaction_roles(gem, medium_table = medium_table)
  validated <- rc_validate_gem(gem)
  full <- gem
  full$S <- validated$S
  full$lb <- validated$lb
  full$ub <- validated$ub
  if (!is.null(gem$reaction_meta)) {
    full$reaction_meta <- gem$reaction_meta[
      match(
        validated$reactions,
        as.character(gem$reaction_meta$reaction_id)
      ),
      , drop = FALSE
    ]
  }

  reference <- full
  medium_diagnostics <- data.frame()
  medium_contract <- list(
    changed_reactions = character(),
    n_changed_reactions = 0L,
    n_removed_reactions = 0L
  )
  if (!is.null(medium_table)) {
    applied <- rc_apply_medium_constraints(
      full,
      medium_table,
      condition = condition,
      strict = FALSE
    )
    full <- applied$gem
    medium_diagnostics <- applied$medium_diagnostics
    medium_contract <- .rc_assert_medium_bounds_only(
      reference,
      full,
      medium_table = medium_table,
      context = "COMPASS-style full-GEM medium application"
    )
  }

  constrained <- rc_validate_gem(full)
  if (!identical(validated$reactions, constrained$reactions)) {
    stop(
      "COMPASS-style full-GEM construction must retain every reference ",
      "reaction after medium application.",
      call. = FALSE
    )
  }
  full$reaction_roles <- full$reaction_meta[
    , intersect(
      c(
        "reaction_id", "role", "role_source",
        "role_confidence"
      ),
      colnames(full$reaction_meta)
    ),
    drop = FALSE
  ]
  full$medium_diagnostics <- medium_diagnostics
  full$medium_contract <- medium_contract
  full$flux_consistency_diagnostics <- data.frame()
  full$closure_diagnostics <- data.frame()
  full$target_status <- "not_prechecked"
  full$build_params <- list(
    strategy = "compass_medium_constrained_full_gem",
    context_specific_reconstruction = FALSE,
    reaction_evidence_used = FALSE,
    fastcore_executed = FALSE,
    corda2_executed = FALSE,
    medium_applied = !is.null(medium_table),
    medium_handling = "exchange_bounds_only_no_reaction_deletion",
    medium_direct_reaction_deletion = FALSE,
    target_feasibility =
      "directional_vmax_then_skip_penalty_lp_when_infeasible",
    flux_consistency_pruning = FALSE,
    flux_consistency_algorithm = "none",
    n_input_reactions = length(validated$reactions),
    n_reactions = ncol(full$S),
    n_medium_removed_reactions = 0L,
    n_medium_bound_changes = medium_contract$n_changed_reactions,
    n_metabolites = nrow(full$S)
  )
  full
}

.rc_full_gem_cache_fingerprint <- function(gem) {
  validated <- rc_validate_gem(gem)
  info <- gem$model_info %||% list()
  payload <- list(
    species = as.character(info$species %||% NA_character_),
    source = as.character(info$source %||% NA_character_),
    version = as.character(
      info$model_version %||% info$version %||% NA_character_
    ),
    commit = as.character(
      info$source_commit %||% info$commit %||% NA_character_
    ),
    checksum = as.character(info$checksum %||% NA_character_),
    S = validated$S,
    lb = validated$lb,
    ub = validated$ub
  )
  file <- tempfile("RegCompassR-gem-fingerprint-", fileext = ".rds")
  on.exit(unlink(file, force = TRUE), add = TRUE)
  saveRDS(payload, file, version = 2)
  unname(tools::md5sum(file)[[1L]])
}

.rc_full_gem_medium_fingerprint <- function(medium) {
  payload <- if (is.null(medium)) {
    list(no_constraints = TRUE)
  } else {
    required <- c("exchange_reaction_id", "lb", "ub", "available")
    missing <- setdiff(required, colnames(medium))
    if (length(missing)) {
      stop(
        "Medium fingerprint input is missing: ",
        paste(missing, collapse = ", "),
        ".", call. = FALSE
      )
    }
    columns <- intersect(
      c(
        "exchange_reaction_id", "condition", "available", "lb", "ub",
        ".no_constraints"
      ),
      colnames(medium)
    )
    value <- medium[, columns, drop = FALSE]
    order_columns <- intersect(
      c("condition", "exchange_reaction_id"),
      colnames(value)
    )
    if (length(order_columns) && nrow(value)) {
      value <- value[do.call(order, value[order_columns]), , drop = FALSE]
    }
    rownames(value) <- NULL
    value
  }
  file <- tempfile("RegCompassR-full-gem-medium-", fileext = ".rds")
  on.exit(unlink(file, force = TRUE), add = TRUE)
  saveRDS(payload, file, version = 2)
  unname(tools::md5sum(file)[[1L]])
}

#' Cache one complete medium-constrained full GEM per medium scenario
.rc_build_full_gem_cache_core <- function(
    gem, dirs, medium_scenarios,
    cache_dir = tempfile("RegCompassR_full_gem_cache_"),
    force = FALSE, conditions = NULL,
    solver = NULL, time_limit = NULL,
    flux_consistency_epsilon = NULL) {
  if (!is.data.frame(dirs) ||
      !all(c("reaction_id", "target_direction") %in% colnames(dirs))) {
    stop("`dirs` must contain `reaction_id` and `target_direction`.", call. = FALSE)
  }
  requested_solver <- if (is.null(solver)) NA_character_ else as.character(solver)
  requested_time_limit <- if (is.null(time_limit)) NA_real_ else as.numeric(time_limit)
  requested_epsilon <- if (is.null(flux_consistency_epsilon)) {
    NA_real_
  } else {
    as.numeric(flux_consistency_epsilon)
  }

  medium_scenarios <- .rc_normalize_medium_scenarios(medium_scenarios)
  if (is.null(conditions)) {
    conditions <- if ("condition" %in% colnames(medium_scenarios)) {
      unique(as.character(medium_scenarios$condition))
    } else {
      "all"
    }
    conditions <- setdiff(conditions, c(NA_character_, "", "all"))
    if (!length(conditions)) conditions <- "all"
  }
  conditions <- unique(trimws(as.character(conditions)))
  conditions <- conditions[!is.na(conditions) & nzchar(conditions)]
  if (!length(conditions)) conditions <- "all"

  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  gem_fingerprint <- .rc_full_gem_cache_fingerprint(gem)
  scenarios <- unique(as.character(medium_scenarios$medium_scenario_id))
  combinations <- expand.grid(
    medium_scenario = scenarios,
    condition = conditions,
    stringsAsFactors = FALSE
  )
  model_files <- list()
  summaries <- vector("list", nrow(combinations))
  safe <- function(value) {
    paste(
      sprintf("%02x", as.integer(charToRaw(enc2utf8(value)))),
      collapse = ""
    )
  }

  for (i in seq_len(nrow(combinations))) {
    scenario <- combinations$medium_scenario[[i]]
    condition <- combinations$condition[[i]]
    identity <- paste(scenario, condition, sep = "::")
    medium <- medium_scenarios[
      as.character(medium_scenarios$medium_scenario_id) == scenario,
      , drop = FALSE
    ]
    if (!nrow(medium) ||
        (".no_constraints" %in% colnames(medium) &&
         all(medium$.no_constraints))) {
      medium <- NULL
    }
    medium_fingerprint <- .rc_full_gem_medium_fingerprint(medium)
    file <- file.path(
      cache_dir,
      paste0(
        "full_gem__gem_", gem_fingerprint,
        "__medium_bounds_", medium_fingerprint,
        "__medium_", safe(scenario),
        "__condition_", safe(condition), ".rds"
      )
    )
    rebuild <- !file.exists(file) || isTRUE(force)
    if (!rebuild) {
      full <- tryCatch(readRDS(file), error = function(error) NULL)
      cached_gem_fingerprint <- if (is.list(full)) {
        full$cache_identity$gem_fingerprint %||% NA_character_
      } else {
        NA_character_
      }
      cached_medium_fingerprint <- if (is.list(full)) {
        full$cache_identity$medium_fingerprint %||% NA_character_
      } else {
        NA_character_
      }
      cached_strategy <- if (is.list(full)) {
        full$build_params$strategy %||% NA_character_
      } else {
        NA_character_
      }
      rebuild <- !identical(cached_gem_fingerprint, gem_fingerprint) ||
        !identical(cached_medium_fingerprint, medium_fingerprint) ||
        !identical(cached_strategy, "compass_medium_constrained_full_gem")
    }
    if (rebuild) {
      full <- rc_build_full_gem(
        gem = gem,
        medium_table = medium,
        condition = if (identical(condition, "all")) NULL else condition
      )
      full$condition <- condition
      full$cache_identity <- list(
        gem_fingerprint = gem_fingerprint,
        medium_fingerprint = medium_fingerprint,
        species = gem$model_info$species %||% NA_character_,
        source = gem$model_info$source %||% NA_character_,
        version = gem$model_info$model_version %||%
          gem$model_info$version %||% NA_character_,
        commit = gem$model_info$source_commit %||%
          gem$model_info$commit %||% NA_character_,
        checksum = gem$model_info$checksum %||% NA_character_
      )
      saveRDS(full, file)
    }
    model_files[[identity]] <- file
    build <- full$build_params
    summaries[[i]] <- data.frame(
      cache_key = paste(
        "full_gem", gem_fingerprint, medium_fingerprint,
        scenario, condition, sep = "::"
      ),
      gem_fingerprint = gem_fingerprint,
      medium_fingerprint = medium_fingerprint,
      medium_scenario = scenario,
      condition = condition,
      file = file,
      n_input_reactions = build$n_input_reactions,
      n_reactions = ncol(full$S),
      n_medium_removed_reactions = 0L,
      n_medium_bound_changes = build$n_medium_bound_changes,
      n_metabolites = nrow(full$S),
      medium_applied = isTRUE(build$medium_applied),
      medium_handling = build$medium_handling,
      requested_solver = requested_solver,
      requested_completion_time_limit = requested_time_limit,
      requested_flux_threshold = requested_epsilon,
      build_strategy = "compass_medium_constrained_full_gem",
      target_status = full$target_status,
      model_version = gem$model_info$model_version %||%
        gem$model_info$version %||% NA_character_,
      model_commit = gem$model_info$source_commit %||%
        gem$model_info$commit %||% NA_character_,
      stringsAsFactors = FALSE
    )
  }

  cache <- list()
  for (i in seq_len(nrow(dirs))) {
    reaction <- as.character(dirs$reaction_id[[i]])
    for (j in seq_len(nrow(combinations))) {
      scenario <- combinations$medium_scenario[[j]]
      condition <- combinations$condition[[j]]
      identity <- paste(scenario, condition, sep = "::")
      key <- paste0(
        "reaction=", utils::URLencode(reaction, reserved = TRUE),
        "::direction=", utils::URLencode(
          as.character(dirs$target_direction[[i]]),
          reserved = TRUE
        ),
        "::medium=", utils::URLencode(scenario, reserved = TRUE),
        "::condition=", utils::URLencode(condition, reserved = TRUE)
      )
      cache[[key]] <- list(
        reaction_id = reaction,
        target_direction = as.character(dirs$target_direction[[i]]),
        medium_scenario = scenario,
        condition = condition,
        file = model_files[[identity]],
        build_strategy = "compass_medium_constrained_full_gem"
      )
    }
  }
  attr(cache, "summary") <- do.call(rbind, summaries)
  attr(cache, "structural_scope") <- "medium_x_complete_full_gem"
  attr(cache, "completion_method") <- "none"
  attr(cache, "fastcore_executed") <- FALSE
  attr(cache, "corda2_executed") <- FALSE
  attr(cache, "medium_handling") <-
    "exchange_bounds_only_no_reaction_deletion"
  cache
}

# Progress-aware entry point; the algorithm remains in the core above.
rc_build_full_gem <- function(...) {
  progress_state <- get0(
    ".rc_layer2_progress_state", mode = "environment", inherits = TRUE
  )
  if (!isTRUE(progress_state$active) &&
      is.null(progress_state$current_task)) {
    return(do.call(.rc_build_full_gem_core, list(...)))
  }
  if (!is.null(progress_state$current_task)) {
    return(do.call(.rc_build_full_gem_core, list(...)))
  }
  args <- list(...)
  context <- .rc_layer2_task_context(
    "ALL", .rc_layer2_medium_id(args$medium_table), "full_gem"
  )
  parts_dir <- .rc_layer2_progress_dir_from_cache(
    .rc_layer2_progress_cache_dir_from_frames()
  )
  .rc_layer2_task_event(
    context, "full_gem_medium_bounds", 1L, 2L,
    "applying medium bounds without reaction deletion",
    scope = "structural", run_kind = "primary", parts_dir = parts_dir
  )
  answer <- do.call(.rc_build_full_gem_core, args)
  .rc_layer2_task_event(
    context, "full_gem_model_ready", 2L, 2L,
    detail = paste0("reactions=", ncol(answer$S)),
    scope = "structural", run_kind = "primary", status = "complete",
    parts_dir = parts_dir
  )
  answer
}

# Progress-aware entry point; the algorithm remains in the core above.
rc_build_full_gem_cache <- function(...) {
  progress_state <- get0(
    ".rc_layer2_progress_state", mode = "environment", inherits = TRUE
  )
  answer <- do.call(
    .rc_build_full_gem_cache_core,
    list(...)
  )
  if (identical(progress_state$run_kind, "primary")) {
    summary <- attr(answer, "summary")
    .rc_layer2_overall_event(
      "structural_models_complete", 3L,
      detail = paste0(
        "full_gem_medium_models=",
        if (is.data.frame(summary)) nrow(summary) else 1L
      )
    )
  }
  answer
}
