rc_parse_microcompass_row_id <- function(x) {
  x <- as.character(x)
  required <- c("reaction", "direction", "medium")

  parse_one <- function(id) {
    fields <- strsplit(id, "::", fixed = TRUE)[[1L]]
    equals <- regexpr("=", fields, fixed = TRUE)
    if (any(equals < 2L)) {
      stop(
        "microCOMPASS row IDs must use `reaction=...::direction=...::medium=...`.",
        call. = FALSE
      )
    }
    keys <- substring(fields, 1L, equals - 1L)
    values <- utils::URLdecode(substring(fields, equals + 1L))
    required_counts <- table(factor(keys, levels = required))
    required_values <- values[match(required, keys)]
    invalid <- any(required_counts != 1L) ||
      anyNA(required_values) ||
      any(!nzchar(trimws(required_values))) ||
      !required_values[[2L]] %in% c("forward", "reverse")
    if (invalid) {
      stop(
        paste(
          "microCOMPASS row IDs must contain exactly one non-empty",
          "`reaction`, `direction`, and `medium` field; `direction` must be",
          "`forward` or `reverse`."
        ),
        call. = FALSE
      )
    }
    named <- stats::setNames(values, keys)
    value <- function(name) {
      hit <- named[name]
      if (!length(hit) || is.na(hit[[1L]]) ||
          !nzchar(trimws(hit[[1L]]))) {
        return(NA_character_)
      }
      unname(hit[[1L]])
    }
    data.frame(
      sample_id = value("sample"),
      module_id = value("module"),
      reaction_id = value("reaction"),
      target_direction = value("direction"),
      medium_scenario = value("medium"),
      condition = value("condition"),
      stringsAsFactors = FALSE
    )
  }

  rows <- lapply(x, parse_one)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Prepare direction-specific reaction targets from signed GEM bounds
rc_prepare_directional_targets <- function(
    gem, target_reactions,
    target_direction = c("both", "forward", "reverse"),
    bound_tolerance = 1e-12) {
  target_direction <- match.arg(target_direction)
  if (!is.numeric(bound_tolerance) || length(bound_tolerance) != 1L ||
      !is.finite(bound_tolerance) || bound_tolerance < 0) {
    stop(
      "`bound_tolerance` must be one finite non-negative number.",
      call. = FALSE
    )
  }
  validated <- rc_validate_gem(gem)
  requested <- unique(trimws(as.character(target_reactions)))
  requested <- requested[!is.na(requested) & nzchar(requested)]
  missing <- setdiff(requested, validated$reactions)
  if (length(missing)) {
    stop(
      "Target reactions missing from GEM: ",
      paste(utils::head(missing, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  rows <- lapply(requested, function(reaction) {
    lb <- validated$lb[[reaction]]
    ub <- validated$ub[[reaction]]
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
    directions <- switch(
      target_direction,
      both = c(
        if (forward_allowed) "forward",
        if (reverse_allowed) "reverse"
      ),
      forward = if (forward_allowed) "forward" else character(),
      reverse = if (reverse_allowed) "reverse" else character()
    )
    if (!length(directions)) directions <- "none"
    data.frame(
      reaction_id = reaction,
      target_direction = directions,
      direction_class = direction_class,
      requested_direction = target_direction,
      direction_status = ifelse(
        directions == "none", "no_allowed_direction", "allowed"
      ),
      stringsAsFactors = FALSE
    )
  })
  if (!length(rows)) {
    return(data.frame(
      reaction_id = character(),
      target_direction = character(),
      direction_class = character(),
      requested_direction = character(),
      direction_status = character(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}

.rc_normalize_medium_scenarios <- function(medium_scenarios) {
  if (is.null(medium_scenarios) ||
      (is.data.frame(medium_scenarios) && !nrow(medium_scenarios))) {
    return(data.frame(
      medium_scenario_id = "base",
      exchange_reaction_id = NA_character_,
      lb = NA_real_,
      ub = NA_real_,
      available = FALSE,
      .no_constraints = TRUE,
      stringsAsFactors = FALSE
    ))
  }
  if (!is.data.frame(medium_scenarios)) {
    stop("`medium_scenarios` must be a data.frame.", call. = FALSE)
  }
  if (!"medium_scenario_id" %in% colnames(medium_scenarios)) {
    medium_scenarios$medium_scenario_id <- "custom"
  }
  if (!".no_constraints" %in% colnames(medium_scenarios)) {
    medium_scenarios$.no_constraints <- FALSE
  } else {
    medium_scenarios$.no_constraints <- as.logical(
      medium_scenarios$.no_constraints
    )
    medium_scenarios$.no_constraints[
      is.na(medium_scenarios$.no_constraints)
    ] <- FALSE
  }
  medium_scenarios
}

.rc_validate_shared_medium <- function(medium_scenarios) {
  medium_scenarios <- .rc_normalize_medium_scenarios(medium_scenarios)
  if ("condition" %in% colnames(medium_scenarios)) {
    condition <- trimws(as.character(medium_scenarios$condition))
    condition <- unique(
      condition[!is.na(condition) & nzchar(condition) & condition != "all"]
    )
    if (length(condition)) {
      stop(
        paste(
          "Shared-GEM scoring requires condition-invariant medium constraints;",
          "remove condition-specific medium rows."
        ),
        call. = FALSE
      )
    }
  }
  medium_scenarios
}

.rc_cache_gem <- function(entry) {
  if (is.list(entry) && !is.null(entry$file)) readRDS(entry$file) else entry
}
