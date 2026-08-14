.rc_resolved_medium_bounds_fingerprint <- function(reference_gem, constrained_gem) {
  reference <- rc_validate_gem(reference_gem)
  constrained <- rc_validate_gem(constrained_gem)
  if (!identical(reference$reactions, constrained$reactions) ||
      !identical(reference$metabolites, constrained$metabolites) ||
      !identical(reference$S, constrained$S)) {
    stop(
      "Resolved-medium fingerprinting requires identical GEM structure; only bounds may differ.",
      call. = FALSE
    )
  }
  payload <- list(
    schema_version = "regcompass_resolved_medium_bounds_v1",
    base_gem_fingerprint = .rc_full_gem_cache_fingerprint(reference_gem),
    bounds = data.frame(
      reaction_id = as.character(reference$reactions),
      final_lb = as.numeric(constrained$lb),
      final_ub = as.numeric(constrained$ub),
      stringsAsFactors = FALSE
    )
  )
  file <- tempfile("RegCompassR-resolved-medium-", fileext = ".rds")
  on.exit(unlink(file, force = TRUE), add = TRUE)
  saveRDS(payload, file, version = 2)
  unname(tools::md5sum(file)[[1L]])
}

.rc_resolved_medium_bound_diagnostics <- function(reference_gem, constrained_gem) {
  reference <- rc_validate_gem(reference_gem)
  constrained <- rc_validate_gem(constrained_gem)
  if (!identical(reference$reactions, constrained$reactions) ||
      !identical(reference$metabolites, constrained$metabolites) ||
      !identical(reference$S, constrained$S)) {
    stop(
      "Resolved-medium diagnostics require identical GEM structure.",
      call. = FALSE
    )
  }
  changed <- reference$lb != constrained$lb | reference$ub != constrained$ub
  list(
    n_changed_bounds_vs_reference = sum(changed),
    resolved_bounds_identical_to_reference = !any(changed),
    changed_reactions = reference$reactions[changed]
  )
}

# Replace the input-table cache identity with the actual LP feasible-region
# identity. The public rc_build_full_gem_cache() wrapper resolves this symbol at
# runtime, so collating this file after compass_medium_semantics.R is sufficient.
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

    # Resolve the actual GEM bounds before deriving cache identity. Medium
    # application is only a deterministic bound transformation here; no CORDA2
    # or LP solve is performed by this resolution step.
    resolved <- rc_build_full_gem(
      gem = gem,
      medium_table = medium,
      condition = if (identical(condition, "all")) NULL else condition
    )
    resolved$condition <- condition
    medium_fingerprint <- .rc_resolved_medium_bounds_fingerprint(gem, resolved)
    bound_diagnostics <- .rc_resolved_medium_bound_diagnostics(gem, resolved)
    input_fingerprint <- .rc_compass_medium_fingerprint(medium)

    file <- file.path(
      cache_dir,
      paste0(
        "full_gem__gem_", gem_fingerprint,
        "__resolved_bounds_", medium_fingerprint,
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
        full$cache_identity$resolved_medium_fingerprint %||%
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
      full <- resolved
      full$cache_identity <- list(
        gem_fingerprint = gem_fingerprint,
        medium_fingerprint = medium_fingerprint,
        resolved_medium_fingerprint = medium_fingerprint,
        medium_input_fingerprint = input_fingerprint,
        medium_fingerprint_semantics =
          "base_GEM_identity_plus_canonical_final_reaction_bounds",
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
    resolved <- NULL
    model_files[[identity]] <- file
    build <- full$build_params
    summaries[[i]] <- data.frame(
      cache_key = paste(
        "full_gem", gem_fingerprint, medium_fingerprint,
        scenario, condition, sep = "::"
      ),
      gem_fingerprint = gem_fingerprint,
      medium_fingerprint = medium_fingerprint,
      resolved_medium_fingerprint = medium_fingerprint,
      medium_input_fingerprint = input_fingerprint,
      medium_fingerprint_semantics =
        "base_GEM_identity_plus_canonical_final_reaction_bounds",
      medium_scenario = scenario,
      condition = condition,
      file = file,
      n_input_reactions = build$n_input_reactions,
      n_reactions = ncol(full$S),
      n_medium_removed_reactions = 0L,
      n_medium_bound_changes = build$n_medium_bound_changes,
      n_changed_bounds_vs_reference =
        bound_diagnostics$n_changed_bounds_vs_reference,
      resolved_bounds_identical_to_reference =
        bound_diagnostics$resolved_bounds_identical_to_reference,
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
          as.character(dirs$target_direction[[i]]), reserved = TRUE
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
  attr(cache, "medium_fingerprint_semantics") <-
    "base_GEM_identity_plus_canonical_final_reaction_bounds"
  cache
}
