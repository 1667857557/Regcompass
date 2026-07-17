# Preserve explicit 0/0 exchange bounds when applying COMPASS-style media.
# `available = FALSE` in the legacy medium application path means "use the
# default secretion bound", not "keep this exchange blocked". Every row in the
# COMPASS-style table therefore remains available and carries its exact capped
# model bounds, including 0/0 blocked exchanges.

.rc_compass_model_bound_medium <- function(gem, exchange_limit = 1) {
  if (!is.numeric(exchange_limit) || length(exchange_limit) != 1L ||
      !is.finite(exchange_limit) || exchange_limit <= 0) {
    stop("`exchange_limit` must be one positive finite number.",
         call. = FALSE)
  }

  validated <- rc_validate_gem(gem)
  if (is.null(gem$reaction_meta) ||
      !"role" %in% colnames(gem$reaction_meta)) {
    gem <- rc_annotate_reaction_roles(gem)
  }
  meta <- gem$reaction_meta[
    match(validated$reactions, as.character(gem$reaction_meta$reaction_id)),
    ,
    drop = FALSE
  ]
  exchange <- as.character(meta$reaction_id[
    as.character(meta$role) == "exchange"
  ])
  exchange <- intersect(exchange, validated$reactions)
  if (!length(exchange)) {
    stop(
      "No `exchange` reactions found for COMPASS-style model bounds.",
      call. = FALSE
    )
  }

  index <- match(exchange, validated$reactions)
  original_lb <- as.numeric(validated$lb[index])
  original_ub <- as.numeric(validated$ub[index])
  lb <- pmax(original_lb, -exchange_limit)
  ub <- pmin(original_ub, exchange_limit)
  if (any(lb > ub)) {
    stop("Capping exchange bounds produced invalid lower/upper bounds.",
         call. = FALSE)
  }

  data.frame(
    medium_scenario_id = "compass_model_bounds",
    exchange_reaction_id = exchange,
    metabolite_id = if ("metabolite_id" %in% colnames(meta)) {
      as.character(meta$metabolite_id[index])
    } else {
      NA_character_
    },
    condition = "all",
    lb = lb,
    ub = ub,
    available = TRUE,
    original_lb = original_lb,
    original_ub = original_ub,
    exchange_limit = exchange_limit,
    evidence_source =
      "gem_directionality_with_compass_uniform_exchange_limit",
    assumption_level = "shared_model_defined_environment",
    target_exchange_flag = FALSE,
    concentration_used_for_rate_bound = FALSE,
    rate_bound_source = "gem_bounds_capped_like_compass",
    stringsAsFactors = FALSE
  )
}
