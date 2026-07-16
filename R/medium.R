#' Build condition-aware medium scenarios
#' @export
rc_make_medium_scenarios <- function(
    gem,
    scenario = c(
      "blood_like", "minimal", "culture_like", "tumor_low_glucose",
      "low_glucose", "low_glutamine", "lactate_available", "custom"
    ),
    custom_medium = NULL,
    uptake_scale = c(1, 0.5, 0.1),
    condition_col = NULL,
    exchange_roles = c("exchange"),
    condition = condition_col) {
  scenario <- match.arg(scenario, several.ok = TRUE)
  if (!is.null(condition) &&
      (length(condition) != 1L || is.na(condition) || !nzchar(as.character(condition)))) {
    stop("`condition` must be NULL or one non-empty condition value.", call. = FALSE)
  }
  if ("custom" %in% scenario) {
    if (is.null(custom_medium)) {
      stop("`custom_medium` is required when `scenario` includes 'custom'.", call. = FALSE)
    }
    required <- c(
      "medium_scenario_id", "exchange_reaction_id", "lb", "ub", "available"
    )
    missing <- setdiff(required, colnames(custom_medium))
    if (length(missing)) {
      stop(
        "`custom_medium` missing columns: ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
  }

  validated <- rc_validate_gem(gem)
  if (is.null(gem$reaction_meta) || !"role" %in% colnames(gem$reaction_meta)) {
    gem <- rc_annotate_reaction_roles(gem)
  }
  meta <- gem$reaction_meta[
    match(validated$reactions, as.character(gem$reaction_meta$reaction_id)),
    , drop = FALSE
  ]
  exchange_roles <- unique(trimws(as.character(exchange_roles)))
  exchange_roles <- exchange_roles[!is.na(exchange_roles) & nzchar(exchange_roles)]
  if (!length(exchange_roles)) {
    stop("`exchange_roles` must contain at least one role.", call. = FALSE)
  }
  exchange <- as.character(meta$reaction_id[
    as.character(meta$role) %in% exchange_roles
  ])
  exchange <- intersect(exchange, validated$reactions)
  built_in <- setdiff(scenario, "custom")
  if (!length(exchange) && length(built_in)) {
    n_boundary <- sum(as.character(meta$role) == "boundary_like", na.rm = TRUE)
    stop(
      "No `exchange` reactions found in GEM reaction metadata. ",
      "`rc_make_medium_scenarios()` requires curated or reliably inferred exchange reactions. ",
      "Found ", n_boundary, " `boundary_like` reactions.",
      call. = FALSE
    )
  }

  scale_defaults <- c(
    blood_like = 1,
    culture_like = 1,
    minimal = 0.1,
    tumor_low_glucose = 0.5,
    low_glucose = 0.1,
    low_glutamine = 0.1,
    lactate_available = 1
  )
  if (!is.numeric(uptake_scale) || !length(uptake_scale) ||
      any(!is.finite(uptake_scale)) || any(uptake_scale < 0)) {
    stop("`uptake_scale` must contain finite non-negative values.", call. = FALSE)
  }
  if (!is.null(names(uptake_scale)) && any(nzchar(names(uptake_scale)))) {
    named <- uptake_scale[nzchar(names(uptake_scale))]
    known <- intersect(names(named), names(scale_defaults))
    scale_defaults[known] <- named[known]
  } else if (length(uptake_scale) == 1L) {
    scale_defaults[] <- uptake_scale[[1L]]
  } else if (length(uptake_scale) > 1L) {
    stop(
      paste(
        "`uptake_scale` must be a single global scale or a named vector",
        "using scenario IDs."
      ),
      call. = FALSE
    )
  }

  text_columns <- intersect(
    c(
      "reaction_id", "reaction_name", "name", "description", "equation",
      "metabolite_id", "metabolite_name"
    ),
    colnames(meta)
  )
  annotation_text <- if (length(text_columns)) {
    tolower(apply(meta[, text_columns, drop = FALSE], 1L, paste, collapse = " "))
  } else {
    tolower(as.character(meta$reaction_id))
  }
  names(annotation_text) <- as.character(meta$reaction_id)
  target_pattern <- list(
    tumor_low_glucose = "glucose|d[- ]?glucose|glc",
    low_glucose = "glucose|d[- ]?glucose|glc",
    low_glutamine = "glutamine|gln",
    lactate_available = "lactate|lac"
  )

  make_rows <- function(scenario_id) {
    base_scale <- if (scenario_id == "minimal") {
      scale_defaults[["minimal"]]
    } else {
      scale_defaults[["blood_like"]]
    }
    reaction_scale <- stats::setNames(rep(base_scale, length(exchange)), exchange)
    target_reactions <- character()
    if (scenario_id %in% names(target_pattern)) {
      target_reactions <- exchange[
        grepl(target_pattern[[scenario_id]], annotation_text[exchange], perl = TRUE)
      ]
      if (!length(target_reactions)) {
        stop(
          "Scenario `", scenario_id,
          "` requires a matching annotated exchange reaction. Supply `custom_medium` when model annotations cannot identify it.",
          call. = FALSE
        )
      }
      reaction_scale[target_reactions] <- scale_defaults[[scenario_id]]
    }
    data.frame(
      medium_scenario_id = scenario_id,
      exchange_reaction_id = exchange,
      metabolite_id = if ("metabolite_id" %in% colnames(meta)) {
        as.character(meta$metabolite_id[match(exchange, meta$reaction_id)])
      } else {
        NA_character_
      },
      condition = if (is.null(condition)) "all" else as.character(condition),
      lb = -10 * as.numeric(reaction_scale[exchange]),
      ub = 1000,
      available = TRUE,
      evidence_source = if (scenario_id %in% names(target_pattern)) {
        "annotation_matched_sensitivity_scenario"
      } else {
        "generic_medium_assumption"
      },
      assumption_level = if (scenario_id %in% c("blood_like", "culture_like")) {
        "generic_baseline"
      } else {
        "sensitivity_scenario"
      },
      target_exchange_flag = exchange %in% target_reactions,
      stringsAsFactors = FALSE
    )
  }

  output <- if (length(built_in)) {
    do.call(rbind, lapply(built_in, make_rows))
  } else {
    NULL
  }
  if ("custom" %in% scenario) {
    custom <- custom_medium
    optional <- c(
      "metabolite_id", "condition", "evidence_source", "assumption_level",
      "target_exchange_flag"
    )
    for (name in setdiff(optional, colnames(custom))) custom[[name]] <- NA
    if (is.null(output)) {
      output <- custom
    } else {
      missing_in_custom <- setdiff(colnames(output), colnames(custom))
      for (name in missing_in_custom) custom[[name]] <- NA
      output <- rbind(output, custom[, colnames(output), drop = FALSE])
    }
  }
  rownames(output) <- NULL
  output
}

#' Apply a medium table to exchange-reaction bounds
rc_apply_medium_constraints <- function(
    gem, medium_table, condition = NULL,
    exchange_default_lb = 0, exchange_default_ub = 1000,
    allow_secretion = TRUE, strict = TRUE) {
  if (!is.logical(allow_secretion) || length(allow_secretion) != 1L ||
      is.na(allow_secretion) || !is.logical(strict) || length(strict) != 1L ||
      is.na(strict)) {
    stop("`allow_secretion` and `strict` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.finite(exchange_default_lb) || !is.finite(exchange_default_ub) ||
      exchange_default_lb > exchange_default_ub) {
    stop("Default exchange bounds must be finite and ordered.", call. = FALSE)
  }
  validated <- rc_validate_gem(gem)
  if (is.null(gem$reaction_meta) || !"role" %in% colnames(gem$reaction_meta)) {
    gem <- rc_annotate_reaction_roles(gem, medium_table = medium_table)
  }
  meta <- gem$reaction_meta[
    match(validated$reactions, as.character(gem$reaction_meta$reaction_id)),
    , drop = FALSE
  ]
  is_exchange <- as.character(meta$role) == "exchange"
  old_lb <- validated$lb
  old_ub <- validated$ub
  lb <- old_lb
  ub <- old_ub
  lb[is_exchange] <- exchange_default_lb
  ub[is_exchange] <- if (allow_secretion) {
    exchange_default_ub
  } else {
    min(exchange_default_ub, 0)
  }
  status <- rep("not_exchange", length(lb))
  status[is_exchange] <- "exchange_default_closed"

  if (!is.null(medium_table)) {
    if (!is.data.frame(medium_table)) {
      stop("`medium_table` must be a data.frame.", call. = FALSE)
    }
    required <- c("exchange_reaction_id", "lb", "ub", "available")
    missing <- setdiff(required, colnames(medium_table))
    if (length(missing)) {
      stop(
        "`medium_table` missing columns: ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
    medium <- medium_table
    medium$exchange_reaction_id <- trimws(as.character(medium$exchange_reaction_id))
    medium$condition <- if ("condition" %in% colnames(medium)) {
      as.character(medium$condition)
    } else {
      "all"
    }
    medium$condition[is.na(medium$condition) | !nzchar(medium$condition)] <- "all"
    keep <- medium$condition == "all" |
      (!is.null(condition) & medium$condition == as.character(condition))
    medium <- medium[keep, , drop = FALSE]
    if (nrow(medium)) {
      if (anyNA(medium$exchange_reaction_id) ||
          any(!nzchar(medium$exchange_reaction_id))) {
        stop("Medium exchange reaction IDs must be non-empty.", call. = FALSE)
      }
      medium$available <- as.logical(medium$available)
      medium$lb <- suppressWarnings(as.numeric(medium$lb))
      medium$ub <- suppressWarnings(as.numeric(medium$ub))
      if (anyNA(medium$available) || any(!is.finite(medium$lb)) ||
          any(!is.finite(medium$ub)) || any(medium$lb > medium$ub)) {
        stop("Medium rows require logical availability and finite ordered bounds.", call. = FALSE)
      }
      priority <- ifelse(medium$condition == "all", 1L, 2L)
      order_index <- order(priority, seq_len(nrow(medium)))
      medium <- medium[order_index, , drop = FALSE]
      duplicate_key <- paste(medium$exchange_reaction_id, medium$condition, sep = "\001")
      if (anyDuplicated(duplicate_key)) {
        stop(
          "`medium_table` contains duplicated reaction/condition rows.",
          call. = FALSE
        )
      }
      unknown <- setdiff(medium$exchange_reaction_id, validated$reactions)
      if (length(unknown)) {
        message <- paste(
          "Medium exchange reactions missing from GEM:",
          paste(utils::head(unknown, 10L), collapse = ", ")
        )
        if (strict) stop(message, call. = FALSE) else warning(message, call. = FALSE)
      }
      medium <- medium[
        medium$exchange_reaction_id %in% validated$reactions,
        , drop = FALSE
      ]
      non_exchange <- medium$exchange_reaction_id[
        !is_exchange[match(medium$exchange_reaction_id, validated$reactions)]
      ]
      if (length(non_exchange)) {
        message <- paste(
          "Medium rows reference reactions not annotated as exchange:",
          paste(utils::head(unique(non_exchange), 10L), collapse = ", ")
        )
        if (strict) stop(message, call. = FALSE) else warning(message, call. = FALSE)
      }
      for (i in seq_len(nrow(medium))) {
        reaction <- medium$exchange_reaction_id[[i]]
        if (!isTRUE(medium$available[[i]])) {
          lb[[reaction]] <- exchange_default_lb
          ub[[reaction]] <- if (allow_secretion) exchange_default_ub else min(exchange_default_ub, 0)
          status[[reaction]] <- "medium_unavailable"
          next
        }
        lb[[reaction]] <- medium$lb[[i]]
        ub[[reaction]] <- if (allow_secretion) {
          medium$ub[[i]]
        } else {
          min(medium$ub[[i]], 0)
        }
        status[[reaction]] <- "medium_available"
      }
    }
  }
  if (any(lb > ub)) {
    stop("Applied medium constraints produced lower bounds above upper bounds.", call. = FALSE)
  }
  gem$lb <- lb
  gem$ub <- ub
  gem$medium_policy <- "condition_aware_exchange_bounds"
  diagnostics <- data.frame(
    reaction_id = names(lb),
    old_lb = old_lb,
    old_ub = old_ub,
    new_lb = lb,
    new_ub = ub,
    medium_status = status,
    condition = condition %||% "all",
    stringsAsFactors = FALSE
  )
  list(gem = gem, medium_diagnostics = diagnostics)
}
