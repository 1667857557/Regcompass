#' Apply condition-aware medium constraints to exchange reactions
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
#' @export
rc_layer2_apply_bounds <- function(gem, medium_policy = NULL) { .Deprecated("rc_apply_medium_constraints"); rc_apply_medium_constraints(gem, medium_policy, strict = FALSE)$gem }
