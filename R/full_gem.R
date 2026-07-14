#' Build a COMPASS-style full GEM for one medium scenario
#'
#' The complete validated Human-GEM is retained. Optional medium constraints are
#' applied once and the resulting model is reused for every target and unit.
#' @export
rc_build_full_gem <- function(gem, medium_table = NULL, condition = NULL) {
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
  full$closure_diagnostics <- data.frame()
  full$target_status <- "not_prechecked"
  full$build_params <- list(
    strategy = "full_gem",
    n_reactions = ncol(full$S),
    n_metabolites = nrow(full$S)
  )
  full
}

#' Cache one complete full GEM per medium scenario
#' @export
rc_build_full_gem_cache <- function(gem, dirs, medium_scenarios,
                                    cache_dir = tempfile(
                                      "RegCompassR_full_gem_cache_"
                                    ),
                                    force = FALSE,
                                    conditions = NULL) {
  if (!is.data.frame(dirs) ||
      !all(c("reaction_id", "target_direction") %in% colnames(dirs))) {
    stop("`dirs` must contain `reaction_id` and `target_direction`.", call. = FALSE)
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
  scenarios <- unique(as.character(medium_scenarios$medium_scenario_id))
  combinations <- expand.grid(
    medium_scenario = scenarios,
    condition = conditions,
    stringsAsFactors = FALSE
  )
  model_files <- list()
  summaries <- vector("list", nrow(combinations))
  safe <- function(value) paste(sprintf("%02x", as.integer(charToRaw(enc2utf8(value)))), collapse = "")

  for (i in seq_len(nrow(combinations))) {
    scenario <- combinations$medium_scenario[[i]]
    condition <- combinations$condition[[i]]
    identity <- paste(scenario, condition, sep = "::")
    file <- file.path(
      cache_dir,
      paste0("full_gem__medium_", safe(scenario), "__condition_", safe(condition), ".rds")
    )
    medium <- medium_scenarios[
      as.character(medium_scenarios$medium_scenario_id) == scenario,
      , drop = FALSE
    ]
    if (!nrow(medium) ||
        (".no_constraints" %in% colnames(medium) && all(medium$.no_constraints))) {
      medium <- NULL
    }
    if (!file.exists(file) || isTRUE(force)) {
      full <- rc_build_full_gem(
        gem = gem,
        medium_table = medium,
        condition = if (identical(condition, "all")) NULL else condition
      )
      full$condition <- condition
      saveRDS(full, file)
    } else {
      full <- readRDS(file)
    }
    model_files[[identity]] <- file
    summaries[[i]] <- data.frame(
      cache_key = paste("full_gem", scenario, condition, sep = "::"),
      medium_scenario = scenario,
      condition = condition,
      file = file,
      n_reactions = ncol(full$S),
      n_metabolites = nrow(full$S),
      build_strategy = "full_gem",
      target_status = "not_prechecked",
      model_version = gem$model_info$model_version %||% gem$model_info$version %||% NA_character_,
      model_commit = gem$model_info$source_commit %||% gem$model_info$commit %||% NA_character_,
      stringsAsFactors = FALSE
    )
  }

  cache <- list()
  for (i in seq_len(nrow(dirs))) {
    for (j in seq_len(nrow(combinations))) {
      scenario <- combinations$medium_scenario[[j]]
      condition <- combinations$condition[[j]]
      identity <- paste(scenario, condition, sep = "::")
      key <- paste(
        dirs$reaction_id[[i]], dirs$target_direction[[i]],
        scenario, condition, sep = "::"
      )
      cache[[key]] <- list(
        reaction_id = as.character(dirs$reaction_id[[i]]),
        target_direction = as.character(dirs$target_direction[[i]]),
        medium_scenario = scenario,
        condition = condition,
        file = model_files[[identity]],
        build_strategy = "full_gem"
      )
    }
  }
  attr(cache, "summary") <- do.call(rbind, summaries)
  cache
}
