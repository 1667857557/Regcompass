# Build one final union GEM per medium scenario.

.rc_safe_cache_token <- function(value) {
  paste(
    sprintf("%02x", as.integer(charToRaw(enc2utf8(as.character(value))))),
    collapse = ""
  )
}

.rc_build_medium_specific_union_gem_cache <- function(
    gem, reaction_membership, core_reactions,
    target_reactions = NULL, medium_scenarios = NULL,
    cache_dir = tempfile("RegCompassR_medium_union_gem_"),
    target_direction = c("both", "forward", "reverse"),
    solver = "highs", time_limit = 300,
    fastcore_epsilon = 1e-4,
    max_support_reactions = 2000,
    strict = TRUE) {
  target_direction <- match.arg(target_direction)
  if (!is.data.frame(reaction_membership) ||
      !"reaction_id" %in% colnames(reaction_membership)) {
    stop("`reaction_membership` must contain reaction_id.", call. = FALSE)
  }
  if (!is.data.frame(core_reactions) ||
      !"reaction_id" %in% colnames(core_reactions)) {
    stop("`core_reactions` must contain reaction_id.", call. = FALSE)
  }
  if ("is_core" %in% colnames(core_reactions)) {
    core_reactions <- core_reactions[
      core_reactions$is_core %in% TRUE, , drop = FALSE
    ]
  }
  if (!is.null(target_reactions)) {
    core_reactions <- core_reactions[
      as.character(core_reactions$reaction_id) %in%
        as.character(target_reactions),
      , drop = FALSE
    ]
  }
  if (!nrow(core_reactions)) {
    stop("No merged core reactions remain for union-GEM scoring.",
         call. = FALSE)
  }

  medium_scenarios <- .rc_validate_shared_medium(medium_scenarios)
  scenarios <- unique(as.character(medium_scenarios$medium_scenario_id))
  scenarios <- scenarios[!is.na(scenarios) & nzchar(scenarios)]
  if (!length(scenarios)) {
    stop("No medium scenarios are available for union-GEM construction.",
         call. = FALSE)
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  cache <- list()
  summaries <- vector("list", length(scenarios))
  names(summaries) <- scenarios
  for (scenario in scenarios) {
    medium <- medium_scenarios[
      as.character(medium_scenarios$medium_scenario_id) == scenario,
      , drop = FALSE
    ]
    if (!nrow(medium) ||
        (".no_constraints" %in% colnames(medium) &&
         all(medium$.no_constraints))) {
      medium <- NULL
    }

    model <- .rc_complete_medium_union_gem(
      gem = gem,
      reaction_membership = reaction_membership,
      core_reactions = core_reactions,
      medium_table = medium,
      target_direction = target_direction,
      solver = solver,
      time_limit = time_limit,
      fastcore_epsilon = fastcore_epsilon,
      max_support_reactions = max_support_reactions,
      strict = strict
    )
    model$shared_across_units <- TRUE
    model$is_union_gem <- TRUE
    model$union_gem_medium_scenario <- scenario
    model$union_gem_scope <-
      "one_medium_shared_across_conditions_and_metacells"

    file <- file.path(
      cache_dir,
      paste0("union_gem__medium_", .rc_safe_cache_token(scenario), ".rds")
    )
    saveRDS(model, file)
    file_checksum <- unname(tools::md5sum(file))
    summaries[[scenario]] <- data.frame(
      medium_scenario = scenario,
      file = file,
      file_checksum = file_checksum,
      n_reactions = ncol(model$S),
      n_metabolites = nrow(model$S),
      n_merged_biological_reactions =
        model$build_params$n_merged_biological_reactions,
      n_global_fastcore_support_reactions =
        model$build_params$n_global_fastcore_support_reactions,
      target_status = model$target_status,
      build_strategy = "medium_specific_union_gem",
      completion_stage =
        "single_global_fastcore_after_meta_module_merge",
      stringsAsFactors = FALSE
    )

    if (!nrow(model$target_directions)) next
    for (i in seq_len(nrow(model$target_directions))) {
      reaction <- as.character(model$target_directions$reaction_id[[i]])
      direction <- as.character(
        model$target_directions$target_direction[[i]]
      )
      key <- paste0(
        "reaction=", utils::URLencode(reaction, reserved = TRUE),
        "::direction=", direction,
        "::medium=", utils::URLencode(scenario, reserved = TRUE)
      )
      cache[[key]] <- list(
        sample_id = "global",
        module_id = "MEDIUM_UNION_GEM",
        reaction_id = reaction,
        target_direction = direction,
        medium_scenario = scenario,
        condition = "all",
        file = file,
        file_checksum = file_checksum,
        build_strategy = "medium_specific_union_gem"
      )
    }
  }
  attr(cache, "summary") <- .rc_bind_frames_fill(summaries)
  cache
}

.rc_read_cached_union_gem <- function(
    file, medium_scenario, expected_checksum = NA_character_) {
  if (!is.character(file) || length(file) != 1L ||
      is.na(file) || !nzchar(file) || !file.exists(file)) {
    stop("A required Stage 5 union GEM cache file is unavailable.",
         call. = FALSE)
  }
  observed_checksum <- unname(tools::md5sum(file))
  if (!is.na(expected_checksum) && nzchar(expected_checksum) &&
      !identical(observed_checksum, as.character(expected_checksum))) {
    stop("A Stage 5 union GEM cache file failed its checksum check.",
         call. = FALSE)
  }
  model <- readRDS(file)
  if (!isTRUE(model$is_union_gem) ||
      !identical(
        as.character(model$union_gem_medium_scenario),
        as.character(medium_scenario)
      )) {
    stop(
      "Second-pass scoring requires the original final medium-specific union GEM.",
      call. = FALSE
    )
  }
  model
}
