# Preserve public validation order, named diagnostics and exact medium-bound
# application after the required result-level corrections.

.rc_required_result_run_microcompass <- rc_run_microcompass
.rc_required_result_apply_medium_constraints <- rc_apply_medium_constraints

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

# The legacy implementation assigns explicit rows with character `[[` indexing.
# Some GEM bound vectors lose or reorder names during conversion, so that
# assignment can append a new element instead of replacing the reaction's bound.
# Reapply validated medium rows by integer reaction index before the GEM enters
# any LP.
rc_apply_medium_constraints <- function(
    gem, medium_table, condition = NULL,
    exchange_default_lb = 0, exchange_default_ub = 1000,
    allow_secretion = TRUE, strict = TRUE) {
  answer <- .rc_required_result_apply_medium_constraints(
    gem = gem,
    medium_table = medium_table,
    condition = condition,
    exchange_default_lb = exchange_default_lb,
    exchange_default_ub = exchange_default_ub,
    allow_secretion = allow_secretion,
    strict = strict
  )
  if (is.null(medium_table) || !is.data.frame(medium_table) ||
      !nrow(medium_table)) {
    return(answer)
  }

  reactions <- colnames(answer$gem$S)
  n_reactions <- length(reactions)
  lb <- as.numeric(answer$gem$lb)[seq_len(n_reactions)]
  ub <- as.numeric(answer$gem$ub)[seq_len(n_reactions)]
  names(lb) <- names(ub) <- reactions

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
    (!is.null(condition) & medium$condition == as.character(condition))
  medium <- medium[keep, , drop = FALSE]

  if (nrow(medium)) {
    medium$available <- as.logical(medium$available)
    medium$lb <- suppressWarnings(as.numeric(medium$lb))
    medium$ub <- suppressWarnings(as.numeric(medium$ub))
    priority <- ifelse(medium$condition == "all", 1L, 2L)
    medium <- medium[order(priority, seq_len(nrow(medium))), , drop = FALSE]

    for (i in seq_len(nrow(medium))) {
      reaction_index <- match(
        medium$exchange_reaction_id[[i]],
        reactions
      )
      if (is.na(reaction_index)) next
      if (!isTRUE(medium$available[[i]])) {
        lb[[reaction_index]] <- exchange_default_lb
        ub[[reaction_index]] <- if (allow_secretion) {
          exchange_default_ub
        } else {
          min(exchange_default_ub, 0)
        }
      } else {
        lb[[reaction_index]] <- medium$lb[[i]]
        ub[[reaction_index]] <- if (allow_secretion) {
          medium$ub[[i]]
        } else {
          min(medium$ub[[i]], 0)
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
  answer$gem$lb <- lb
  answer$gem$ub <- ub
  if (is.data.frame(answer$medium_diagnostics)) {
    index <- match(
      answer$medium_diagnostics$reaction_id,
      reactions
    )
    answer$medium_diagnostics$new_lb <- lb[index]
    answer$medium_diagnostics$new_ub <- ub[index]
  }
  answer
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
