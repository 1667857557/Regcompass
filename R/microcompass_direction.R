#' Prepare direction-specific reaction targets from signed Human-GEM bounds
#' @export
rc_prepare_directional_targets <- function(gem, target_reactions,
                                           target_direction = c("both", "forward", "reverse"),
                                           bound_tolerance = 1e-12) {
  target_direction <- match.arg(target_direction)
  gv <- rc_validate_gem(gem)
  target_reactions <- intersect(unique(as.character(target_reactions)), gv$reactions)
  if (!length(target_reactions)) {
    return(data.frame(
      reaction_id = character(), target_direction = character(),
      direction_class = character(), stringsAsFactors = FALSE
    ))
  }

  rows <- lapply(target_reactions, function(r) {
    lb <- gv$lb[[r]]
    ub <- gv$ub[[r]]
    forward_allowed <- ub > bound_tolerance
    reverse_allowed <- lb < -bound_tolerance
    direction_class <- if (forward_allowed && reverse_allowed) {
      "reversible"
    } else if (forward_allowed) {
      "forward_only"
    } else if (reverse_allowed) {
      "reverse_only"
    } else {
      "blocked"
    }
    dirs <- switch(
      target_direction,
      both = c(if (forward_allowed) "forward", if (reverse_allowed) "reverse"),
      forward = if (forward_allowed) "forward" else character(),
      reverse = if (reverse_allowed) "reverse" else character()
    )
    if (!length(dirs)) return(NULL)
    data.frame(
      reaction_id = r,
      target_direction = dirs,
      direction_class = direction_class,
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) {
    return(data.frame(
      reaction_id = character(), target_direction = character(),
      direction_class = character(), stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}

.rc_normalize_medium_scenarios <- function(medium_scenarios) {
  if (is.null(medium_scenarios) || (is.data.frame(medium_scenarios) && !nrow(medium_scenarios))) {
    return(data.frame(
      medium_scenario_id = "base",
      exchange_reaction_id = NA_character_,
      lb = NA_real_, ub = NA_real_, available = FALSE,
      .no_constraints = TRUE,
      stringsAsFactors = FALSE
    ))
  }
  if (!is.data.frame(medium_scenarios)) stop("`medium_scenarios` must be a data.frame.", call. = FALSE)
  if (!"medium_scenario_id" %in% colnames(medium_scenarios)) {
    medium_scenarios$medium_scenario_id <- "custom"
  }
  medium_scenarios$.no_constraints <- FALSE
  medium_scenarios
}

.rc_cache_gem <- function(entry) {
  if (is.list(entry) && !is.null(entry$file)) readRDS(entry$file) else entry
}

.rc_unit_sample_map <- function(unit_meta, unit_ids, sample_col) {
  if (is.null(unit_meta) || !sample_col %in% colnames(unit_meta)) {
    stop("Meta-module GEM mode requires unit metadata with `", sample_col, "`.", call. = FALSE)
  }
  id_cols <- intersect(c("unit_id", "pool_id", "metacell_id"), colnames(unit_meta))
  for (id_col in id_cols) {
    out <- stats::setNames(as.character(unit_meta[[sample_col]]), as.character(unit_meta[[id_col]]))
    if (all(unit_ids %in% names(out))) return(out[unit_ids])
  }
  if (nrow(unit_meta) == length(unit_ids)) {
    return(stats::setNames(as.character(unit_meta[[sample_col]]), unit_ids))
  }
  stop("Could not align Layer 2 units to sample IDs.", call. = FALSE)
}
