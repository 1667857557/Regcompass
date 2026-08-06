# COMPASS-compatible extracellular exchange-bound semantics -----------------
#
# COMPASS first applies a uniform cap to the uptake direction of every
# exchange reaction, preserves the secretion direction, and then lets explicit
# medium entries override the uptake cap. A reaction omitted from the medium
# therefore retains the capped model-defined uptake direction instead of being
# closed. This file is collated after the legacy medium implementation and
# installs the corrected semantics without changing the public medium table or
# Layer 2 interfaces.

.rc_compass_medium_semantics_version <- "compass_exchange_bounds_v2"

.rc_exchange_uptake_directions <- function(gem, reaction_ids, strict = TRUE) {
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("`strict` must be TRUE or FALSE.", call. = FALSE)
  }
  validated <- rc_validate_gem(gem)
  ids <- unique(trimws(as.character(reaction_ids)))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  unknown <- setdiff(ids, validated$reactions)
  if (length(unknown)) {
    message <- paste(
      "Exchange reactions missing from GEM:",
      paste(utils::head(unknown, 10L), collapse = ", ")
    )
    if (strict) stop(message, call. = FALSE) else warning(message, call. = FALSE)
  }
  ids <- intersect(ids, validated$reactions)
  direction <- vapply(ids, function(id) {
    index <- match(id, validated$reactions)
    coefficients <- as.numeric(validated$S[, index, drop = FALSE])
    coefficients <- coefficients[is.finite(coefficients) & coefficients != 0]
    if (!length(coefficients)) return("none")
    if (all(coefficients < 0)) return("reverse")
    if (all(coefficients > 0)) return("forward")
    "mixed"
  }, character(1))
  unresolved <- direction %in% c("none", "mixed")
  if (any(unresolved) && strict) {
    stop(
      paste0(
        "Exchange uptake direction is not one-sided for: ",
        paste(utils::head(ids[unresolved], 10L), collapse = ", "),
        ". COMPASS exchange limits require all non-zero stoichiometric ",
        "coefficients to have the same sign."
      ),
      call. = FALSE
    )
  }
  data.frame(
    reaction_id = ids,
    uptake_direction = unname(direction),
    stringsAsFactors = FALSE
  )
}

.rc_compass_validate_limit <- function(value, name, allow_infinite = FALSE) {
  value <- suppressWarnings(as.numeric(value))
  valid <- length(value) == 1L && !is.na(value) && value >= 0 &&
    (is.finite(value) || (allow_infinite && is.infinite(value) && value > 0))
  if (!valid) {
    stop(
      "`", name, "` must be one non-negative ",
      if (allow_infinite) "numeric value." else "finite numeric value.",
      call. = FALSE
    )
  }
  value
}

.rc_compass_baseline_bounds <- function(
    old_lb, old_ub, uptake_direction, exchange_limit,
    unlisted_policy = c("compass", "closed"),
    allow_secretion = TRUE, secretion_limit = Inf) {
  unlisted_policy <- match.arg(unlisted_policy)
  exchange_limit <- .rc_compass_validate_limit(exchange_limit, "exchange_limit")
  secretion_limit <- .rc_compass_validate_limit(
    secretion_limit, "exchange_default_ub", allow_infinite = TRUE
  )
  lb <- old_lb
  ub <- old_ub
  reverse <- names(uptake_direction)[uptake_direction == "reverse"]
  forward <- names(uptake_direction)[uptake_direction == "forward"]
  unresolved <- names(uptake_direction)[
    !uptake_direction %in% c("reverse", "forward")
  ]

  if (identical(unlisted_policy, "compass")) {
    lb[reverse] <- pmax(old_lb[reverse], -exchange_limit)
    ub[forward] <- pmin(old_ub[forward], exchange_limit)
  } else {
    lb[reverse] <- pmax(old_lb[reverse], 0)
    ub[forward] <- pmin(old_ub[forward], 0)
    if (length(unresolved)) {
      lb[unresolved] <- pmax(old_lb[unresolved], 0)
      ub[unresolved] <- pmin(old_ub[unresolved], 0)
    }
  }

  if (isTRUE(allow_secretion)) {
    if (is.finite(secretion_limit)) {
      ub[reverse] <- pmin(ub[reverse], secretion_limit)
      lb[forward] <- pmax(lb[forward], -secretion_limit)
    }
  } else {
    ub[reverse] <- pmin(ub[reverse], 0)
    lb[forward] <- pmax(lb[forward], 0)
    if (length(unresolved)) {
      lb[unresolved] <- pmax(lb[unresolved], 0)
      ub[unresolved] <- pmin(ub[unresolved], 0)
    }
  }
  list(lb = lb, ub = ub)
}

.rc_compass_medium_row_scope <- function(medium, index) {
  if ("bound_scope" %in% colnames(medium)) {
    scope <- tolower(trimws(as.character(medium$bound_scope[[index]])))
    if (!is.na(scope) && nzchar(scope)) {
      if (!scope %in% c("uptake", "both")) {
        stop("Medium `bound_scope` must be `uptake` or `both`.", call. = FALSE)
      }
      return(scope)
    }
  }
  evidence <- if ("evidence_source" %in% colnames(medium)) {
    as.character(medium$evidence_source[[index]])
  } else {
    NA_character_
  }
  assumption <- if ("assumption_level" %in% colnames(medium)) {
    as.character(medium$assumption_level[[index]])
  } else {
    NA_character_
  }
  if ((!is.na(evidence) && evidence %in% c(
    "literature_backed_medium_catalog",
    "gem_directionality_with_uniform_exchange_cap"
  )) || (!is.na(assumption) && assumption %in% c(
    "availability_catalog_with_relative_uptake_cap",
    "shared_model_defined_environment"
  ))) {
    return("uptake")
  }
  "both"
}

.rc_compass_medium_row_uptake_limit <- function(
    medium, index, uptake_direction) {
  if ("uptake_limit" %in% colnames(medium)) {
    value <- suppressWarnings(as.numeric(medium$uptake_limit[[index]]))
    if (length(value) == 1L && is.finite(value) && value >= 0) return(value)
  }
  if (all(c("exchange_limit", "uptake_fraction") %in% colnames(medium))) {
    cap <- suppressWarnings(as.numeric(medium$exchange_limit[[index]]))
    fraction <- suppressWarnings(as.numeric(medium$uptake_fraction[[index]]))
    if (is.finite(cap) && cap >= 0 && is.finite(fraction) && fraction >= 0) {
      return(cap * fraction)
    }
  }
  if ("exchange_limit" %in% colnames(medium)) {
    cap <- suppressWarnings(as.numeric(medium$exchange_limit[[index]]))
    evidence <- if ("evidence_source" %in% colnames(medium)) {
      as.character(medium$evidence_source[[index]])
    } else {
      NA_character_
    }
    if (is.finite(cap) && cap >= 0 &&
        !is.na(evidence) &&
        identical(evidence, "gem_directionality_with_uniform_exchange_cap")) {
      return(cap)
    }
  }
  if (identical(uptake_direction, "reverse")) {
    value <- -suppressWarnings(as.numeric(medium$lb[[index]]))
  } else if (identical(uptake_direction, "forward")) {
    value <- suppressWarnings(as.numeric(medium$ub[[index]]))
  } else {
    value <- NA_real_
  }
  if (length(value) == 1L && is.finite(value) && value >= 0) value else NA_real_
}

.rc_compass_close_uptake <- function(
    old_lb, old_ub, uptake_direction,
    allow_secretion = TRUE, secretion_limit = Inf) {
  lb <- old_lb
  ub <- old_ub
  if (identical(uptake_direction, "reverse")) {
    lb <- max(old_lb, 0)
    if (!isTRUE(allow_secretion)) {
      ub <- min(old_ub, 0)
    } else if (is.finite(secretion_limit)) {
      ub <- min(old_ub, secretion_limit)
    }
  } else if (identical(uptake_direction, "forward")) {
    ub <- min(old_ub, 0)
    if (!isTRUE(allow_secretion)) {
      lb <- max(old_lb, 0)
    } else if (is.finite(secretion_limit)) {
      lb <- max(old_lb, -secretion_limit)
    }
  } else {
    lb <- max(old_lb, 0)
    ub <- min(old_ub, 0)
  }
  c(lb = lb, ub = ub)
}

.rc_apply_compass_medium_constraints <- function(
    gem, medium_table, condition = NULL,
    exchange_default_lb = NULL, exchange_default_ub = Inf,
    allow_secretion = TRUE, strict = TRUE,
    exchange_limit = 1,
    unlisted_policy = c("compass", "closed")) {
  policy_missing <- missing(unlisted_policy)
  unlisted_policy <- match.arg(unlisted_policy)
  if (!is.logical(allow_secretion) || length(allow_secretion) != 1L ||
      is.na(allow_secretion) ||
      !is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("`allow_secretion` and `strict` must be TRUE or FALSE.", call. = FALSE)
  }
  exchange_limit <- .rc_compass_validate_limit(exchange_limit, "exchange_limit")
  secretion_limit <- .rc_compass_validate_limit(
    exchange_default_ub, "exchange_default_ub", allow_infinite = TRUE
  )
  if (!is.null(exchange_default_lb)) {
    legacy_lb <- suppressWarnings(as.numeric(exchange_default_lb))
    if (length(legacy_lb) != 1L || !is.finite(legacy_lb) || legacy_lb > 0) {
      stop(
        "Legacy `exchange_default_lb` must be zero or a finite negative cap.",
        call. = FALSE
      )
    }
    if (legacy_lb == 0 && policy_missing) {
      unlisted_policy <- "closed"
    } else if (legacy_lb < 0) {
      exchange_limit <- abs(legacy_lb)
    }
  }

  validated <- rc_validate_gem(gem)
  reactions <- validated$reactions
  if (is.null(gem$reaction_meta) || !"role" %in% colnames(gem$reaction_meta)) {
    gem <- rc_annotate_reaction_roles(gem, medium_table = medium_table)
  }
  meta <- gem$reaction_meta[
    match(reactions, as.character(gem$reaction_meta$reaction_id)),
    , drop = FALSE
  ]
  is_exchange <- as.character(meta$role) == "exchange"
  is_exchange[is.na(is_exchange)] <- FALSE
  exchange_ids <- reactions[is_exchange]
  orientation <- .rc_exchange_uptake_directions(
    gem, exchange_ids, strict = strict
  )
  uptake_direction <- stats::setNames(
    rep("not_exchange", length(reactions)), reactions
  )
  uptake_direction[orientation$reaction_id] <- orientation$uptake_direction

  old_lb <- stats::setNames(as.numeric(validated$lb), reactions)
  old_ub <- stats::setNames(as.numeric(validated$ub), reactions)
  exchange_direction <- uptake_direction[exchange_ids]
  baseline <- .rc_compass_baseline_bounds(
    old_lb = old_lb[exchange_ids],
    old_ub = old_ub[exchange_ids],
    uptake_direction = exchange_direction,
    exchange_limit = exchange_limit,
    unlisted_policy = unlisted_policy,
    allow_secretion = allow_secretion,
    secretion_limit = secretion_limit
  )
  lb <- old_lb
  ub <- old_ub
  lb[exchange_ids] <- baseline$lb
  ub[exchange_ids] <- baseline$ub
  baseline_lb <- lb
  baseline_ub <- ub

  status <- stats::setNames(rep("not_exchange", length(reactions)), reactions)
  if (identical(unlisted_policy, "compass")) {
    status[exchange_ids[exchange_direction == "reverse"]] <-
      "compass_unlisted_reverse_uptake_capped"
    status[exchange_ids[exchange_direction == "forward"]] <-
      "compass_unlisted_forward_uptake_capped"
  } else {
    status[exchange_ids] <- "closed_world_unlisted_uptake_closed"
  }
  unresolved_ids <- exchange_ids[
    !exchange_direction %in% c("reverse", "forward")
  ]
  status[unresolved_ids] <- "exchange_uptake_direction_unresolved"
  listed <- stats::setNames(rep(FALSE, length(reactions)), reactions)
  applied_scope <- stats::setNames(rep(NA_character_, length(reactions)), reactions)

  medium <- medium_table
  if (!is.null(medium)) {
    if (!is.data.frame(medium)) {
      stop("`medium_table` must be a data.frame.", call. = FALSE)
    }
    required <- c("exchange_reaction_id", "lb", "ub", "available")
    missing_columns <- setdiff(required, colnames(medium))
    if (length(missing_columns)) {
      stop(
        "`medium_table` missing columns: ",
        paste(missing_columns, collapse = ", "),
        call. = FALSE
      )
    }
    medium$exchange_reaction_id <-
      trimws(as.character(medium$exchange_reaction_id))
    medium$condition <- if ("condition" %in% colnames(medium)) {
      as.character(medium$condition)
    } else {
      "all"
    }
    medium$condition[is.na(medium$condition) | !nzchar(medium$condition)] <-
      "all"
    keep <- medium$condition == "all"
    if (!is.null(condition)) {
      keep <- keep | medium$condition == as.character(condition)
    }
    medium <- medium[keep, , drop = FALSE]

    if (nrow(medium) && "unlisted_policy" %in% colnames(medium)) {
      requested_policy <- unique(tolower(trimws(as.character(
        medium$unlisted_policy
      ))))
      requested_policy <- requested_policy[
        !is.na(requested_policy) & nzchar(requested_policy)
      ]
      if (length(requested_policy) > 1L ||
          (length(requested_policy) &&
           !requested_policy %in% c("compass", "closed"))) {
        stop(
          "Medium `unlisted_policy` must contain at most one of ",
          "`compass` or `closed`.", call. = FALSE
        )
      }
      if (length(requested_policy) &&
          !identical(requested_policy, unlisted_policy)) {
        unlisted_policy <- requested_policy
        baseline <- .rc_compass_baseline_bounds(
          old_lb = old_lb[exchange_ids],
          old_ub = old_ub[exchange_ids],
          uptake_direction = exchange_direction,
          exchange_limit = exchange_limit,
          unlisted_policy = unlisted_policy,
          allow_secretion = allow_secretion,
          secretion_limit = secretion_limit
        )
        lb[exchange_ids] <- baseline$lb
        ub[exchange_ids] <- baseline$ub
        baseline_lb <- lb
        baseline_ub <- ub
        if (identical(unlisted_policy, "compass")) {
          status[exchange_ids[exchange_direction == "reverse"]] <-
            "compass_unlisted_reverse_uptake_capped"
          status[exchange_ids[exchange_direction == "forward"]] <-
            "compass_unlisted_forward_uptake_capped"
        } else {
          status[exchange_ids] <- "closed_world_unlisted_uptake_closed"
        }
        status[unresolved_ids] <- "exchange_uptake_direction_unresolved"
      }
    }

    if (nrow(medium)) {
      medium$available <- as.logical(medium$available)
      medium$lb <- suppressWarnings(as.numeric(medium$lb))
      medium$ub <- suppressWarnings(as.numeric(medium$ub))
      if (anyNA(medium$exchange_reaction_id) ||
          any(!nzchar(medium$exchange_reaction_id)) ||
          anyNA(medium$available) ||
          any(!is.finite(medium$lb)) ||
          any(!is.finite(medium$ub)) ||
          any(medium$lb > medium$ub)) {
        stop(
          "Medium rows require valid reaction IDs, logical availability, ",
          "and finite ordered bounds.", call. = FALSE
        )
      }
      duplicate_key <- paste(
        medium$exchange_reaction_id, medium$condition, sep = "\001"
      )
      if (anyDuplicated(duplicate_key)) {
        stop("`medium_table` contains duplicated reaction/condition rows.",
             call. = FALSE)
      }
      unknown <- setdiff(medium$exchange_reaction_id, reactions)
      if (length(unknown)) {
        message <- paste(
          "Medium exchange reactions missing from GEM:",
          paste(utils::head(unknown, 10L), collapse = ", ")
        )
        if (strict) stop(message, call. = FALSE) else warning(message, call. = FALSE)
      }
      medium <- medium[
        medium$exchange_reaction_id %in% reactions,
        , drop = FALSE
      ]
      reaction_index <- match(medium$exchange_reaction_id, reactions)
      non_exchange <- medium$exchange_reaction_id[!is_exchange[reaction_index]]
      if (length(non_exchange)) {
        message <- paste(
          "Medium rows reference reactions not annotated as exchange:",
          paste(utils::head(unique(non_exchange), 10L), collapse = ", ")
        )
        if (strict) stop(message, call. = FALSE) else warning(message, call. = FALSE)
        keep_exchange <- is_exchange[reaction_index]
        medium <- medium[keep_exchange, , drop = FALSE]
        reaction_index <- reaction_index[keep_exchange]
      }

      for (i in seq_len(nrow(medium))) {
        index <- reaction_index[[i]]
        reaction <- reactions[[index]]
        direction <- uptake_direction[[reaction]]
        scope <- .rc_compass_medium_row_scope(medium, i)
        listed[[reaction]] <- TRUE
        applied_scope[[reaction]] <- scope

        if (!isTRUE(medium$available[[i]])) {
          closed <- .rc_compass_close_uptake(
            old_lb[[reaction]], old_ub[[reaction]], direction,
            allow_secretion = allow_secretion,
            secretion_limit = secretion_limit
          )
          lb[[reaction]] <- closed[["lb"]]
          ub[[reaction]] <- closed[["ub"]]
          status[[reaction]] <- "medium_unavailable_uptake_closed"
          next
        }

        if (identical(scope, "uptake")) {
          uptake_cap <- .rc_compass_medium_row_uptake_limit(
            medium, i, direction
          )
          if (!is.finite(uptake_cap) || uptake_cap < 0) {
            message <- paste0(
              "Cannot resolve an uptake cap for medium reaction `", reaction,
              "`. Supply `uptake_limit`, `exchange_limit` plus ",
              "`uptake_fraction`, or an orientation-compatible bound."
            )
            if (strict) stop(message, call. = FALSE) else warning(message, call. = FALSE)
            scope <- "both"
            applied_scope[[reaction]] <- scope
          } else if (identical(direction, "reverse")) {
            lb[[reaction]] <- max(old_lb[[reaction]], -uptake_cap)
            ub[[reaction]] <- old_ub[[reaction]]
          } else if (identical(direction, "forward")) {
            lb[[reaction]] <- old_lb[[reaction]]
            ub[[reaction]] <- min(old_ub[[reaction]], uptake_cap)
          } else {
            scope <- "both"
            applied_scope[[reaction]] <- scope
          }
          if (identical(scope, "uptake")) {
            if (!isTRUE(allow_secretion)) {
              if (identical(direction, "reverse")) {
                ub[[reaction]] <- min(ub[[reaction]], 0)
              } else {
                lb[[reaction]] <- max(lb[[reaction]], 0)
              }
            } else if (is.finite(secretion_limit)) {
              if (identical(direction, "reverse")) {
                ub[[reaction]] <- min(ub[[reaction]], secretion_limit)
              } else {
                lb[[reaction]] <- max(lb[[reaction]], -secretion_limit)
              }
            }
            status[[reaction]] <- "medium_available_uptake_override"
          }
        }

        if (identical(scope, "both")) {
          lb[[reaction]] <- max(old_lb[[reaction]], medium$lb[[i]])
          ub[[reaction]] <- min(old_ub[[reaction]], medium$ub[[i]])
          if (!isTRUE(allow_secretion)) {
            if (identical(direction, "reverse")) {
              ub[[reaction]] <- min(ub[[reaction]], 0)
            } else if (identical(direction, "forward")) {
              lb[[reaction]] <- max(lb[[reaction]], 0)
            }
          }
          status[[reaction]] <- "medium_available_bounds_intersection"
        }
      }
    }
  }

  if (any(lb > ub)) {
    bad <- reactions[lb > ub]
    stop(
      "Applied medium constraints produced lower bounds above upper bounds for: ",
      paste(utils::head(bad, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  gem$S <- validated$S
  gem$lb <- stats::setNames(as.numeric(lb), reactions)
  gem$ub <- stats::setNames(as.numeric(ub), reactions)
  gem$medium_policy <-
    "compass_uptake_cap_then_explicit_directional_overrides"
  gem$medium_semantics_version <- .rc_compass_medium_semantics_version
  gem$medium_unlisted_policy <- unlisted_policy

  diagnostics <- data.frame(
    reaction_id = reactions,
    old_lb = as.numeric(old_lb),
    old_ub = as.numeric(old_ub),
    baseline_lb = as.numeric(baseline_lb),
    baseline_ub = as.numeric(baseline_ub),
    new_lb = as.numeric(gem$lb),
    new_ub = as.numeric(gem$ub),
    lower_bound_expanded = gem$lb < old_lb,
    upper_bound_expanded = gem$ub > old_ub,
    expanded_relative_to_compass_baseline =
      gem$lb < baseline_lb | gem$ub > baseline_ub,
    uptake_direction = as.character(uptake_direction),
    medium_listed = as.logical(listed),
    bound_scope = as.character(applied_scope),
    medium_status = as.character(status),
    unlisted_policy = unlisted_policy,
    condition = condition %||% "all",
    stringsAsFactors = FALSE
  )
  if (any(diagnostics$lower_bound_expanded | diagnostics$upper_bound_expanded)) {
    stop("Medium application expanded the original GEM feasible region.",
         call. = FALSE)
  }
  list(gem = gem, medium_diagnostics = diagnostics)
}

.rc_make_compass_model_bound_medium <- function(gem, exchange_limit = 1) {
  validated <- rc_validate_gem(gem)
  exchange_limit <- .rc_compass_validate_limit(exchange_limit, "exchange_limit")
  if (is.null(gem$reaction_meta) || !"role" %in% colnames(gem$reaction_meta)) {
    gem <- rc_annotate_reaction_roles(gem)
  }
  meta <- gem$reaction_meta[
    match(validated$reactions, as.character(gem$reaction_meta$reaction_id)),
    , drop = FALSE
  ]
  is_exchange <- as.character(meta$role) == "exchange"
  is_exchange[is.na(is_exchange)] <- FALSE
  exchange <- validated$reactions[is_exchange]
  if (!length(exchange)) {
    stop("No exchange reactions were identified in the GEM.", call. = FALSE)
  }
  orientation <- .rc_exchange_uptake_directions(gem, exchange, strict = TRUE)
  direction <- stats::setNames(
    orientation$uptake_direction, orientation$reaction_id
  )[exchange]
  index <- match(exchange, validated$reactions)
  original_lb <- stats::setNames(as.numeric(validated$lb[index]), exchange)
  original_ub <- stats::setNames(as.numeric(validated$ub[index]), exchange)
  bounds <- .rc_compass_baseline_bounds(
    original_lb, original_ub, direction,
    exchange_limit = exchange_limit,
    unlisted_policy = "compass",
    allow_secretion = TRUE,
    secretion_limit = Inf
  )
  exchange_meta <- meta[
    match(exchange, as.character(meta$reaction_id)), , drop = FALSE
  ]
  exchange_map <- .rc_medium_exchange_metabolites(
    gem, exchange_meta, validated
  )
  data.frame(
    medium_scenario_id = "compass_model_bounds",
    exchange_reaction_id = exchange,
    metabolite_id = exchange_map$metabolite_id,
    gem_metabolite_name = exchange_map$gem_metabolite_name,
    match_method = exchange_map$mapping_source,
    preset_metabolite = NA_character_,
    nutrient_category = NA_character_,
    concentration_mM = NA_real_,
    concentration_basis = NA_character_,
    component_reference_doi = NA_character_,
    condition = "all",
    lb = as.numeric(bounds$lb[exchange]),
    ub = as.numeric(bounds$ub[exchange]),
    available = TRUE,
    original_lb = as.numeric(original_lb[exchange]),
    original_ub = as.numeric(original_ub[exchange]),
    exchange_limit = exchange_limit,
    uptake_fraction = 1,
    uptake_limit = exchange_limit,
    uptake_direction = as.character(direction[exchange]),
    bound_scope = "uptake",
    unlisted_policy = "compass",
    evidence_source = "gem_directionality_with_uniform_exchange_cap",
    assumption_level = "shared_model_defined_environment",
    target_exchange_flag = FALSE,
    concentration_used_for_rate_bound = FALSE,
    rate_bound_source =
      "original_gem_uptake_direction_capped_without_secretion_cap",
    stringsAsFactors = FALSE
  )
}

.rc_assert_compass_medium_bounds_only <- function(
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

  role_gem <- reference_gem
  if (is.null(role_gem$reaction_meta) ||
      !all(c("reaction_id", "role") %in% colnames(role_gem$reaction_meta))) {
    role_gem <- rc_annotate_reaction_roles(
      role_gem, medium_table = medium_table
    )
  }
  meta <- role_gem$reaction_meta[
    match(reference$reactions, as.character(role_gem$reaction_meta$reaction_id)),
    , drop = FALSE
  ]
  is_exchange <- as.character(meta$role) == "exchange"
  is_exchange[is.na(is_exchange)] <- FALSE
  exchange <- reference$reactions[is_exchange]
  changed <- reference$reactions[
    reference$lb != constrained$lb | reference$ub != constrained$ub
  ]
  unexpected <- setdiff(changed, exchange)
  if (length(unexpected)) {
    stop(
      context,
      " modified bounds outside annotated exchange reactions: ",
      paste(utils::head(unexpected, 10L), collapse = ", "),
      ".", call. = FALSE
    )
  }
  supplied <- if (!is.null(medium_table) &&
                  "exchange_reaction_id" %in% colnames(medium_table)) {
    unique(as.character(medium_table$exchange_reaction_id))
  } else {
    character()
  }
  supplied <- supplied[!is.na(supplied) & nzchar(supplied)]
  invisible(list(
    changed_reactions = changed,
    changed_supplied_reactions = intersect(changed, supplied),
    changed_unlisted_exchange_reactions = setdiff(changed, supplied),
    n_changed_reactions = length(changed),
    n_changed_supplied_reactions = length(intersect(changed, supplied)),
    n_changed_unlisted_exchange_reactions = length(setdiff(changed, supplied)),
    n_removed_reactions = 0L,
    medium_semantics_version = .rc_compass_medium_semantics_version
  ))
}

.rc_compass_medium_fingerprint <- function(medium) {
  medium_payload <- if (is.null(medium)) {
    list(no_constraints = TRUE)
  } else {
    required <- c("exchange_reaction_id", "lb", "ub", "available")
    missing_columns <- setdiff(required, colnames(medium))
    if (length(missing_columns)) {
      stop(
        "Medium fingerprint input is missing: ",
        paste(missing_columns, collapse = ", "),
        ".", call. = FALSE
      )
    }
    columns <- intersect(
      c(
        "exchange_reaction_id", "condition", "available", "lb", "ub",
        "bound_scope", "uptake_limit", "uptake_fraction", "exchange_limit",
        "unlisted_policy", ".no_constraints"
      ),
      colnames(medium)
    )
    value <- medium[, columns, drop = FALSE]
    order_columns <- intersect(
      c("condition", "exchange_reaction_id"), colnames(value)
    )
    if (length(order_columns) && nrow(value)) {
      value <- value[do.call(order, value[order_columns]), , drop = FALSE]
    }
    rownames(value) <- NULL
    value
  }
  payload <- list(
    medium_semantics_version = .rc_compass_medium_semantics_version,
    unlisted_default = "compass",
    medium = medium_payload
  )
  file <- tempfile("RegCompassR-compass-medium-", fileext = ".rds")
  on.exit(unlink(file, force = TRUE), add = TRUE)
  saveRDS(payload, file, version = 2)
  unname(tools::md5sum(file)[[1L]])
}

# Install the corrected implementations after the legacy medium files are
# collated. Existing callers keep the same function names and data contracts.
rc_apply_medium_constraints <- .rc_apply_compass_medium_constraints
.rc_compass_model_bound_medium <- .rc_make_compass_model_bound_medium
.rc_assert_medium_bounds_only <- .rc_assert_compass_medium_bounds_only
.rc_full_gem_medium_fingerprint <- .rc_compass_medium_fingerprint
