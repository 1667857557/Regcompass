#' Apply condition-aware medium constraints to exchange reactions
#' @export
rc_make_medium_scenarios <- function(gem,
                                     scenario = c("blood_like", "minimal", "culture_like", "tumor_low_glucose", "low_glucose", "low_glutamine", "lactate_available", "custom"),
                                     custom_medium = NULL,
                                     uptake_scale = c(1, 0.5, 0.1),
                                     condition_col = NULL) {
  scenario <- match.arg(scenario, several.ok = TRUE)
  if ("custom" %in% scenario) {
    if (is.null(custom_medium)) stop("`custom_medium` is required when `scenario` includes 'custom'.", call. = FALSE)
    req <- c("medium_scenario_id", "exchange_reaction_id", "lb", "ub", "available")
    miss <- setdiff(req, colnames(custom_medium)); if (length(miss)) stop("`custom_medium` missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  gv <- rc_validate_gem(gem)
  if (is.null(gem$reaction_meta) || !"role" %in% colnames(gem$reaction_meta)) gem <- rc_annotate_reaction_roles(gem)
  meta <- gem$reaction_meta
  ex <- as.character(meta$reaction_id[as.character(meta$role) == "exchange"])
  ex <- intersect(ex, gv$reactions)
  make_rows <- function(sc) {
    scale <- switch(sc, minimal = 0.1, blood_like = 1, culture_like = 1,
                    tumor_low_glucose = 0.5, low_glucose = 0.5,
                    low_glutamine = 0.5, lactate_available = 1, 1)
    data.frame(
      medium_scenario_id = sc,
      exchange_reaction_id = ex,
      metabolite_id = if ("metabolite_id" %in% colnames(meta)) as.character(meta$metabolite_id[match(ex, meta$reaction_id)]) else NA_character_,
      condition = if (is.null(condition_col)) "all" else as.character(condition_col),
      lb = -10 * scale,
      ub = 1000,
      available = TRUE,
      evidence_source = "curated_scenario_assumption",
      assumption_level = if (sc %in% c("blood_like", "minimal")) "generic_physiological_scenario" else "sensitivity_scenario",
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, lapply(setdiff(scenario, "custom"), make_rows))
  if ("custom" %in% scenario) {
    cm <- custom_medium
    for (nm in setdiff(c("metabolite_id","condition","evidence_source","assumption_level"), colnames(cm))) cm[[nm]] <- NA
    out <- rbind(out, cm[, colnames(out), drop = FALSE])
  }
  rownames(out) <- NULL
  out
}

#' @export
rc_apply_medium_constraints <- function(gem, medium_table, condition = NULL, exchange_default_lb = 0,
                                        exchange_default_ub = 1000, allow_secretion = TRUE, strict = TRUE) {
  gv <- rc_validate_gem(gem)
  if (is.null(gem$reaction_meta) || !"role" %in% colnames(gem$reaction_meta)) gem <- rc_annotate_reaction_roles(gem)
  meta <- gem$reaction_meta[match(gv$reactions, as.character(gem$reaction_meta$reaction_id)), , drop = FALSE]
  is_ex <- as.character(meta$role) == "exchange"
  old_lb <- gv$lb; old_ub <- gv$ub; lb <- old_lb; ub <- old_ub
  lb[is_ex] <- exchange_default_lb
  ub[is_ex] <- if (allow_secretion) pmax(exchange_default_ub, 0) else pmin(exchange_default_ub, 0)
  status <- rep("not_exchange", length(lb)); status[is_ex] <- "exchange_default_closed"
  if (!is.null(medium_table)) {
    req <- c("exchange_reaction_id", "lb", "ub", "available")
    miss <- setdiff(req, colnames(medium_table)); if (length(miss)) stop("`medium_table` missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
    mt <- medium_table
    mt$condition <- if ("condition" %in% colnames(mt)) as.character(mt$condition) else "all"
    keep <- mt$condition == "all" | (!is.null(condition) & mt$condition == condition)
    mt <- mt[keep, , drop = FALSE]
    mt$.prio <- ifelse(mt$condition == "all", 1L, 2L); mt <- mt[order(mt$.prio), , drop = FALSE]
    unk <- setdiff(as.character(mt$exchange_reaction_id), gv$reactions)
    if (length(unk)) { msg <- paste("Medium exchange reactions missing from GEM:", paste(utils::head(unk, 10), collapse = ", ")); if (strict) stop(msg, call. = FALSE) else warning(msg, call. = FALSE) }
    mt <- mt[as.character(mt$exchange_reaction_id) %in% gv$reactions, , drop = FALSE]
    for (i in seq_len(nrow(mt))) {
      r <- as.character(mt$exchange_reaction_id[i]); if (!isTRUE(as.logical(mt$available[i]))) next
      lb[r] <- as.numeric(mt$lb[i]); ub[r] <- as.numeric(mt$ub[i]); status[r] <- "medium_available"
    }
  }
  gem$lb <- lb; gem$ub <- ub; gem$medium_policy <- "condition_aware_exchange_bounds"
  diag <- data.frame(reaction_id = names(lb), old_lb = old_lb, old_ub = old_ub, new_lb = lb, new_ub = ub, medium_status = status, condition = condition %||% "all", stringsAsFactors = FALSE)
  list(gem = gem, medium_diagnostics = diag)
}
