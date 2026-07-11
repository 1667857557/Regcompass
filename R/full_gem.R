#' Build a COMPASS-style full GEM for one medium scenario
#'
#' Applies optional medium constraints to the complete validated GEM while keeping
#' all model reactions available to the downstream COMPASS LP.
#' @export
rc_build_full_gem <- function(gem, medium_table = NULL, condition = NULL) {
  gem <- rc_annotate_reaction_roles(gem, medium_table = medium_table)
  gv <- rc_validate_gem(gem)
  full <- gem
  full$S <- gv$S
  full$lb <- gv$lb
  full$ub <- gv$ub
  if (!is.null(gem$reaction_meta)) {
    full$reaction_meta <- gem$reaction_meta[match(gv$reactions, as.character(gem$reaction_meta$reaction_id)), , drop = FALSE]
  }
  med_diag <- data.frame()
  if (!is.null(medium_table)) {
    app <- rc_apply_medium_constraints(full, medium_table, condition = condition, strict = FALSE)
    full <- app$gem
    med_diag <- app$medium_diagnostics
  }
  full$reaction_roles <- full$reaction_meta[, intersect(c("reaction_id", "role", "role_source", "role_confidence"), colnames(full$reaction_meta)), drop = FALSE]
  full$medium_diagnostics <- med_diag
  full$closure_diagnostics <- data.frame()
  full$build_params <- list(strategy = "full_gem", n_reactions = ncol(full$S), n_metabolites = nrow(full$S))
  full
}

#' Build full-GEM cache once per medium scenario for COMPASS-style analyses
#' @export
rc_build_full_gem_cache <- function(gem,
                                    dirs,
                                    medium_scenarios,
                                    cache_dir = tempfile("RegCompassR_full_gem_cache_"),
                                    force = FALSE) {
  if (!is.data.frame(dirs) || !all(c("reaction_id", "target_direction") %in% colnames(dirs))) {
    stop("`dirs` must contain `reaction_id` and `target_direction`.", call. = FALSE)
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(medium_scenarios)) {
    medium_scenarios <- data.frame(medium_scenario_id = "base", exchange_reaction_id = character(), lb = numeric(), ub = numeric(), available = logical(), stringsAsFactors = FALSE)
  }
  if (!"medium_scenario_id" %in% colnames(medium_scenarios)) medium_scenarios$medium_scenario_id <- "custom"

  scenarios <- unique(as.character(medium_scenarios$medium_scenario_id))
  scenario_files <- list()
  scenario_summary <- vector("list", length(scenarios))
  for (i in seq_along(scenarios)) {
    sc <- scenarios[[i]]
    rds <- file.path(cache_dir, paste0("full_gem__medium_", gsub("[^A-Za-z0-9_.-]+", "_", sc), ".rds"))
    if (!file.exists(rds) || isTRUE(force)) {
      mt <- medium_scenarios[as.character(medium_scenarios$medium_scenario_id) == sc, , drop = FALSE]
      fg <- rc_build_full_gem(gem = gem, medium_table = mt)
      saveRDS(fg, rds)
    } else {
      fg <- readRDS(rds)
    }
    scenario_files[[sc]] <- rds
    scenario_summary[[i]] <- data.frame(cache_key = paste("full_gem", sc, sep = "::"), medium_scenario = sc, file = rds,
                                        n_reactions = ncol(fg$S), n_metabolites = nrow(fg$S),
                                        build_strategy = "full_gem", fallback_from = NA_character_, target_status = "not_prechecked",
                                        model_version = (gem$model_info$model_version %||% gem$model_info$version %||% NA_character_),
                                        model_commit = (gem$model_info$source_commit %||% gem$model_info$commit %||% NA_character_),
                                        stringsAsFactors = FALSE)
  }

  reaction_cache <- list()
  for (i in seq_len(nrow(dirs))) {
    for (sc in scenarios) {
      row_key <- paste(dirs$reaction_id[[i]], dirs$target_direction[[i]], sc, sep = "::")
      reaction_cache[[row_key]] <- list(reaction_id = dirs$reaction_id[[i]], target_direction = dirs$target_direction[[i]],
                                        medium_scenario = sc, file = scenario_files[[sc]], build_strategy = "full_gem")
    }
  }
  attr(reaction_cache, "summary") <- do.call(rbind, scenario_summary)
  reaction_cache
}
