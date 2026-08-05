#' Build a COMPASS-style full GEM for one medium scenario
#'
#' Optional medium constraints are applied once. In Layer 2 full-GEM mode,
#' reactions that cannot carry non-zero steady-state flux after medium
#' application are removed without using expression evidence, FASTCORE or
#' CORDA2. The resulting model is reused for every target and unit.
rc_build_full_gem <- function(
    gem, medium_table = NULL, condition = NULL,
    prune_flux_inconsistent = FALSE,
    solver = "highs", time_limit = 300,
    flux_consistency_epsilon = 1e-8) {
  if (!is.logical(prune_flux_inconsistent) ||
      length(prune_flux_inconsistent) != 1L ||
      is.na(prune_flux_inconsistent)) {
    stop("`prune_flux_inconsistent` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(time_limit) || length(time_limit) != 1L ||
      !is.finite(time_limit) || time_limit <= 0) {
    stop("`time_limit` must be one positive finite number.", call. = FALSE)
  }
  if (!is.numeric(flux_consistency_epsilon) ||
      length(flux_consistency_epsilon) != 1L ||
      !is.finite(flux_consistency_epsilon) ||
      flux_consistency_epsilon <= 0) {
    stop("`flux_consistency_epsilon` must be positive.", call. = FALSE)
  }

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

  medium_diagnostics <- data.frame()
  if (!is.null(medium_table)) {
    applied <- rc_apply_medium_constraints(
      full,
      medium_table,
      condition = condition,
      strict = FALSE
    )
    full <- applied$gem
    medium_diagnostics <- applied$medium_diagnostics
  }

  medium_validated <- rc_validate_gem(full)
  input_reactions <- medium_validated$reactions
  consistent <- input_reactions
  if (isTRUE(prune_flux_inconsistent)) {
    consistent <- .rc_fastcc_consistent_reactions(
      full,
      solver = solver,
      time_limit = time_limit,
      epsilon = flux_consistency_epsilon
    )
    consistent <- input_reactions[input_reactions %in% consistent]
    if (!length(consistent)) {
      stop(
        "The medium-constrained full GEM contains no flux-consistent reactions.",
        call. = FALSE
      )
    }
    full <- .rc_subset_gem(full, consistent)
  }
  removed <- setdiff(input_reactions, consistent)
  flux_consistency_diagnostics <- data.frame(
    reaction_id = input_reactions,
    flux_consistent = input_reactions %in% consistent,
    structural_action = ifelse(
      input_reactions %in% consistent,
      "retained",
      "removed_medium_flux_inconsistent"
    ),
    stringsAsFactors = FALSE
  )

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
  full$flux_consistency_diagnostics <- flux_consistency_diagnostics
  full$closure_diagnostics <- flux_consistency_diagnostics
  full$target_status <- if (length(removed)) {
    "medium_flux_inconsistent_reactions_removed"
  } else {
    "all_reactions_medium_flux_consistent"
  }
  full$build_params <- list(
    strategy = if (isTRUE(prune_flux_inconsistent)) {
      "medium_flux_consistency_pruned_full_gem"
    } else {
      "full_gem"
    },
    context_specific_reconstruction = FALSE,
    reaction_evidence_used = FALSE,
    fastcore_executed = FALSE,
    corda2_executed = FALSE,
    medium_applied = !is.null(medium_table),
    flux_consistency_pruning = isTRUE(prune_flux_inconsistent),
    flux_consistency_algorithm = if (isTRUE(prune_flux_inconsistent)) {
      "FASTCC_flux_consistency_only"
    } else {
      "none"
    },
    flux_consistency_epsilon = flux_consistency_epsilon,
    solver = solver,
    time_limit = time_limit,
    n_input_reactions = length(input_reactions),
    n_reactions = ncol(full$S),
    n_flux_inconsistent_reactions = length(removed),
    removed_flux_inconsistent_reactions = removed,
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

.rc_full_gem_pruning_fingerprint <- function(
    flux_consistency_epsilon, solver) {
  payload <- list(
    algorithm = "medium_flux_consistency_pruned_full_gem_v1",
    epsilon = as.numeric(flux_consistency_epsilon),
    solver = as.character(solver)
  )
  file <- tempfile("RegCompassR-full-gem-pruning-", fileext = ".rds")
  on.exit(unlink(file, force = TRUE), add = TRUE)
  saveRDS(payload, file, version = 2)
  unname(tools::md5sum(file)[[1L]])
}

#' Cache one medium-pruned full GEM per medium scenario
rc_build_full_gem_cache <- function(
    gem, dirs, medium_scenarios,
    cache_dir = tempfile("RegCompassR_full_gem_cache_"),
    force = FALSE, conditions = NULL,
    solver = NULL, time_limit = NULL,
    flux_consistency_epsilon = NULL) {
  if (!is.data.frame(dirs) ||
      !all(c("reaction_id", "target_direction") %in% colnames(dirs))) {
    stop("`dirs` must contain `reaction_id` and `target_direction`.", call. = FALSE)
  }
  context <- .rc_layer2_completion_context
  solver <- as.character(
    solver %||% context$solver %||% "highs"
  )
  time_limit <- as.numeric(
    time_limit %||% context$completion_time_limit %||% 300
  )
  flux_consistency_epsilon <- as.numeric(
    flux_consistency_epsilon %||% context$flux_threshold %||% 1e-8
  )
  if (length(solver) != 1L || is.na(solver) || !nzchar(solver)) {
    stop("`solver` must be one non-empty solver name.", call. = FALSE)
  }
  if (length(time_limit) != 1L || !is.finite(time_limit) ||
      time_limit <= 0) {
    stop("`time_limit` must be one positive finite number.", call. = FALSE)
  }
  if (length(flux_consistency_epsilon) != 1L ||
      !is.finite(flux_consistency_epsilon) ||
      flux_consistency_epsilon <= 0) {
    stop("`flux_consistency_epsilon` must be positive.", call. = FALSE)
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
  pruning_fingerprint <- .rc_full_gem_pruning_fingerprint(
    flux_consistency_epsilon, solver
  )
  scenarios <- unique(as.character(medium_scenarios$medium_scenario_id))
  combinations <- expand.grid(
    medium_scenario = scenarios,
    condition = conditions,
    stringsAsFactors = FALSE
  )
  model_files <- list()
  model_reactions <- list()
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
    file <- file.path(
      cache_dir,
      paste0(
        "full_gem__gem_", gem_fingerprint,
        "__pruning_", pruning_fingerprint,
        "__medium_", safe(scenario),
        "__condition_", safe(condition), ".rds"
      )
    )
    medium <- medium_scenarios[
      as.character(medium_scenarios$medium_scenario_id) == scenario,
      , drop = FALSE
    ]
    if (!nrow(medium) ||
        (".no_constraints" %in% colnames(medium) &&
         all(medium$.no_constraints))) {
      medium <- NULL
    }
    rebuild <- !file.exists(file) || isTRUE(force)
    if (!rebuild) {
      full <- tryCatch(readRDS(file), error = function(error) NULL)
      cached_gem_fingerprint <- if (is.list(full)) {
        full$cache_identity$gem_fingerprint %||% NA_character_
      } else {
        NA_character_
      }
      cached_pruning_fingerprint <- if (is.list(full)) {
        full$cache_identity$pruning_fingerprint %||% NA_character_
      } else {
        NA_character_
      }
      cached_strategy <- if (is.list(full)) {
        full$build_params$strategy %||% NA_character_
      } else {
        NA_character_
      }
      rebuild <- !identical(cached_gem_fingerprint, gem_fingerprint) ||
        !identical(cached_pruning_fingerprint, pruning_fingerprint) ||
        !identical(
          cached_strategy,
          "medium_flux_consistency_pruned_full_gem"
        )
    }
    if (rebuild) {
      full <- rc_build_full_gem(
        gem = gem,
        medium_table = medium,
        condition = if (identical(condition, "all")) NULL else condition,
        prune_flux_inconsistent = TRUE,
        solver = solver,
        time_limit = time_limit,
        flux_consistency_epsilon = flux_consistency_epsilon
      )
      full$condition <- condition
      full$cache_identity <- list(
        gem_fingerprint = gem_fingerprint,
        pruning_fingerprint = pruning_fingerprint,
        solver = solver,
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
    model_reactions[[identity]] <- colnames(full$S)
    build <- full$build_params
    summaries[[i]] <- data.frame(
      cache_key = paste(
        "full_gem", gem_fingerprint, pruning_fingerprint,
        scenario, condition, sep = "::"
      ),
      gem_fingerprint = gem_fingerprint,
      pruning_fingerprint = pruning_fingerprint,
      medium_scenario = scenario,
      condition = condition,
      file = file,
      n_input_reactions = build$n_input_reactions,
      n_reactions = ncol(full$S),
      n_flux_inconsistent_reactions =
        build$n_flux_inconsistent_reactions,
      n_metabolites = nrow(full$S),
      flux_consistency_epsilon = flux_consistency_epsilon,
      solver = solver,
      completion_time_limit = time_limit,
      build_strategy = "medium_flux_consistency_pruned_full_gem",
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
      if (!reaction %in% model_reactions[[identity]]) next
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
        build_strategy = "medium_flux_consistency_pruned_full_gem"
      )
    }
  }
  if (!length(cache)) {
    stop(
      "No requested target reactions remain after medium flux-consistency pruning.",
      call. = FALSE
    )
  }
  attr(cache, "summary") <- do.call(rbind, summaries)
  attr(cache, "structural_scope") <- "medium_x_full_gem"
  attr(cache, "completion_method") <- "none"
  attr(cache, "fastcore_executed") <- FALSE
  attr(cache, "corda2_executed") <- FALSE
  cache
}
