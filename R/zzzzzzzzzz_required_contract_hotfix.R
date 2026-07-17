# Preserve public validation order, named diagnostics and exact medium-bound
# application after the required result-level corrections.

.rc_required_result_run_microcompass <- rc_run_microcompass

rc_run_microcompass <- function(
    layer1, gem, target_reactions = NULL,
    medium_table = NULL, medium_scenarios = NULL,
    mode = c("full_gem", "meta_module_gem"),
    reaction_membership = NULL, core_reactions = NULL,
    unit = c("sample_celltype", "metacell"),
    condition_col = "condition", sample_col = "sample_id",
    celltype_col = "cell_type", model_params = list(),
    penalty_weights = c(expr = 1.0, confidence = 0.5, missing = 1.0),
    omega = 0.95,
    target_direction = c("both", "forward", "reverse"),
    parallel = TRUE,
    solver = c("highs", "gurobi", "glpk"),
    time_limit = 60, flux_threshold = 1e-8,
    BPPARAM = NULL) {
  # Invalid public arguments must fail before attempting to inspect the GEM.
  mode <- match.arg(mode)
  unit <- match.arg(unit)
  target_direction <- match.arg(target_direction)
  solver <- match.arg(solver)

  .rc_required_result_run_microcompass(
    layer1 = layer1,
    gem = gem,
    target_reactions = target_reactions,
    medium_table = medium_table,
    medium_scenarios = medium_scenarios,
    mode = mode,
    reaction_membership = reaction_membership,
    core_reactions = core_reactions,
    unit = unit,
    condition_col = condition_col,
    sample_col = sample_col,
    celltype_col = celltype_col,
    model_params = model_params,
    penalty_weights = penalty_weights,
    omega = omega,
    target_direction = target_direction,
    parallel = parallel,
    solver = solver,
    time_limit = time_limit,
    flux_threshold = flux_threshold,
    BPPARAM = BPPARAM
  )
}

# Bounds are aligned to the validated stoichiometric reaction order and written
# by integer index. This prevents character `[[` assignment from appending a
# duplicate bound entry instead of modifying the reaction used by the LP.
rc_apply_medium_constraints <- function(
    gem, medium_table, condition = NULL,
    exchange_default_lb = 0, exchange_default_ub = 1000,
    allow_secretion = TRUE, strict = TRUE) {
  if (!is.logical(allow_secretion) || length(allow_secretion) != 1L ||
      is.na(allow_secretion) ||
      !is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("`allow_secretion` and `strict` must be TRUE or FALSE.",
         call. = FALSE)
  }
  if (!is.finite(exchange_default_lb) ||
      !is.finite(exchange_default_ub) ||
      exchange_default_lb > exchange_default_ub) {
    stop("Default exchange bounds must be finite and ordered.",
         call. = FALSE)
  }

  validated <- rc_validate_gem(gem)
  reactions <- validated$reactions
  if (is.null(gem$reaction_meta) ||
      !"role" %in% colnames(gem$reaction_meta)) {
    gem <- rc_annotate_reaction_roles(
      gem,
      medium_table = medium_table
    )
  }
  meta <- gem$reaction_meta[
    match(reactions, as.character(gem$reaction_meta$reaction_id)),
    ,
    drop = FALSE
  ]
  is_exchange <- as.character(meta$role) == "exchange"

  old_lb <- stats::setNames(as.numeric(validated$lb), reactions)
  old_ub <- stats::setNames(as.numeric(validated$ub), reactions)
  lb <- old_lb
  ub <- old_ub
  lb[is_exchange] <- exchange_default_lb
  ub[is_exchange] <- if (allow_secretion) {
    exchange_default_ub
  } else {
    min(exchange_default_ub, 0)
  }
  status <- stats::setNames(
    rep("not_exchange", length(reactions)),
    reactions
  )
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
    medium$exchange_reaction_id <- trimws(
      as.character(medium$exchange_reaction_id)
    )
    medium$condition <- if ("condition" %in% colnames(medium)) {
      as.character(medium$condition)
    } else {
      "all"
    }
    medium$condition[
      is.na(medium$condition) | !nzchar(medium$condition)
    ] <- "all"
    keep <- medium$condition == "all" |
      (!is.null(condition) &
         medium$condition == as.character(condition))
    medium <- medium[keep, , drop = FALSE]

    if (nrow(medium)) {
      if (anyNA(medium$exchange_reaction_id) ||
          any(!nzchar(medium$exchange_reaction_id))) {
        stop("Medium exchange reaction IDs must be non-empty.",
             call. = FALSE)
      }
      medium$available <- as.logical(medium$available)
      medium$lb <- suppressWarnings(as.numeric(medium$lb))
      medium$ub <- suppressWarnings(as.numeric(medium$ub))
      if (anyNA(medium$available) ||
          any(!is.finite(medium$lb)) ||
          any(!is.finite(medium$ub)) ||
          any(medium$lb > medium$ub)) {
        stop(
          "Medium rows require logical availability and finite ordered bounds.",
          call. = FALSE
        )
      }

      priority <- ifelse(medium$condition == "all", 1L, 2L)
      medium <- medium[
        order(priority, seq_len(nrow(medium))),
        ,
        drop = FALSE
      ]
      duplicate_key <- paste(
        medium$exchange_reaction_id,
        medium$condition,
        sep = "\001"
      )
      if (anyDuplicated(duplicate_key)) {
        stop(
          "`medium_table` contains duplicated reaction/condition rows.",
          call. = FALSE
        )
      }

      unknown <- setdiff(medium$exchange_reaction_id, reactions)
      if (length(unknown)) {
        message <- paste(
          "Medium exchange reactions missing from GEM:",
          paste(utils::head(unknown, 10L), collapse = ", ")
        )
        if (strict) stop(message, call. = FALSE) else
          warning(message, call. = FALSE)
      }
      medium <- medium[
        medium$exchange_reaction_id %in% reactions,
        ,
        drop = FALSE
      ]

      reaction_index <- match(
        medium$exchange_reaction_id,
        reactions
      )
      non_exchange <- medium$exchange_reaction_id[
        !is_exchange[reaction_index]
      ]
      if (length(non_exchange)) {
        message <- paste(
          "Medium rows reference reactions not annotated as exchange:",
          paste(utils::head(unique(non_exchange), 10L), collapse = ", ")
        )
        if (strict) stop(message, call. = FALSE) else
          warning(message, call. = FALSE)
      }

      for (i in seq_len(nrow(medium))) {
        index <- reaction_index[[i]]
        if (!isTRUE(medium$available[[i]])) {
          lb[[index]] <- exchange_default_lb
          ub[[index]] <- if (allow_secretion) {
            exchange_default_ub
          } else {
            min(exchange_default_ub, 0)
          }
          status[[index]] <- "medium_unavailable"
        } else {
          lb[[index]] <- medium$lb[[i]]
          ub[[index]] <- if (allow_secretion) {
            medium$ub[[i]]
          } else {
            min(medium$ub[[i]], 0)
          }
          status[[index]] <- "medium_available"
        }
      }
    }
  }

  if (any(lb > ub)) {
    stop(
      "Applied medium constraints produced lower bounds above upper bounds.",
      call. = FALSE
    )
  }
  gem$S <- validated$S
  gem$lb <- stats::setNames(as.numeric(lb), reactions)
  gem$ub <- stats::setNames(as.numeric(ub), reactions)
  gem$medium_policy <- "condition_aware_exchange_bounds"

  diagnostics <- data.frame(
    reaction_id = reactions,
    old_lb = as.numeric(old_lb),
    old_ub = as.numeric(old_ub),
    new_lb = as.numeric(gem$lb),
    new_ub = as.numeric(gem$ub),
    medium_status = as.character(status),
    condition = condition %||% "all",
    stringsAsFactors = FALSE
  )
  list(gem = gem, medium_diagnostics = diagnostics)
}

rc_compute_multiome_penalty <- function(...) {
  answer <- .rc_compute_multiome_penalty_core(...)
  reaction_ids <- rownames(answer$penalty)
  named_components <- c(
    "gpr_missing_fraction", "role", "role_source", "role_confidence",
    "role_override_flag", "support_penalty_used"
  )
  for (name in named_components) {
    value <- answer$components[[name]]
    if (!is.null(value) && is.null(dim(value)) &&
        length(value) == length(reaction_ids)) {
      names(value) <- reaction_ids
      answer$components[[name]] <- value
    }
  }
  answer$evidence_policy <- "penalty_only"
  answer$evidence_description <- paste(
    "This is not the original COMPASS expression-neighbourhood penalty.",
    "RegCompass uses a COMPASS-like inverse-support expression term;",
    "multiome evidence modifies the LP objective penalty only and does not",
    "directly change stoichiometry or internal reaction bounds."
  )
  answer
}
