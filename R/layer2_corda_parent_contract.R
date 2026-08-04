# Preserve the original CORDA parent-model semantics without FASTCC pruning.

.rc_fastcore_parent_before_corda_contract <- .rc_fastcore_parent

.rc_corda_parent <- function(
    gem, medium_table = NULL, condition = NULL,
    forbidden_roles = c("demand", "sink", "artificial_support"),
    solver = "highs", time_limit = 300) {
  parent <- rc_build_full_gem(
    gem, medium_table = medium_table, condition = condition
  )
  parent <- rc_annotate_reaction_roles(parent, medium_table = medium_table)
  validated <- rc_validate_gem(parent)
  meta <- parent$reaction_meta[
    match(validated$reactions, as.character(parent$reaction_meta$reaction_id)),
    , drop = FALSE
  ]
  role <- if ("role" %in% colnames(meta)) {
    as.character(meta$role)
  } else {
    rep("unknown", nrow(meta))
  }
  forbidden <- validated$reactions[role %in% forbidden_roles]
  if (length(forbidden)) {
    parent$lb[forbidden] <- 0
    parent$ub[forbidden] <- 0
  }
  validated <- rc_validate_gem(parent)
  feasibility <- rc_solve_lp(
    obj = rep(0, length(validated$reactions)),
    A = validated$S,
    lhs = rep(0, nrow(validated$S)),
    rhs = rep(0, nrow(validated$S)),
    lb = validated$lb,
    ub = validated$ub,
    solver = solver,
    time_limit = time_limit
  )
  if (!identical(feasibility$status, "optimal")) {
    stop(
      "The medium-constrained CORDA parent GEM is not feasible: ",
      feasibility$status,
      call. = FALSE
    )
  }
  parent$corda_forbidden_reactions <- forbidden
  parent$corda_parent_original_lb <- validated$lb
  parent$corda_parent_original_ub <- validated$ub
  parent$corda_parent_prepruning <- "none"
  parent$corda_parent_contract <- paste(
    "medium and RegCompass forbidden-role bounds are applied, but no",
    "FASTCC/FASTCORE consistency pruning is performed before CORDA"
  )
  parent
}

.rc_fastcore_parent <- function(
    gem, medium_table = NULL, condition = NULL,
    forbidden_roles = c("demand", "sink", "artificial_support"),
    solver = "highs", time_limit = 300, fastcore_epsilon = 1e-4) {
  if (isTRUE(getOption("RegCompassR.corda_parent_active", FALSE))) {
    return(.rc_corda_parent(
      gem = gem,
      medium_table = medium_table,
      condition = condition,
      forbidden_roles = forbidden_roles,
      solver = solver,
      time_limit = time_limit
    ))
  }
  .rc_fastcore_parent_before_corda_contract(
    gem = gem,
    medium_table = medium_table,
    condition = condition,
    forbidden_roles = forbidden_roles,
    solver = solver,
    time_limit = time_limit,
    fastcore_epsilon = fastcore_epsilon
  )
}

.rc_complete_celltype_medium_corda_gem_parent_base <-
  .rc_complete_celltype_medium_corda_gem

.rc_complete_celltype_medium_corda_gem <- function(...) {
  args <- list(...)
  previous <- getOption("RegCompassR.corda_parent_active", FALSE)
  options(RegCompassR.corda_parent_active = TRUE)
  on.exit(
    options(RegCompassR.corda_parent_active = previous),
    add = TRUE
  )
  model <- do.call(
    .rc_complete_celltype_medium_corda_gem_parent_base,
    args
  )
  build <- model$build_params
  build$n_parent_reactions <- ncol(model$S) +
    sum(!as.character(model$reaction_meta$reaction_id) %in%
          as.character(model$reaction_meta$reaction_id))
  build$n_fastcc_consistent_parent_reactions <- NULL
  build$n_fastcc_inconsistent_parent_reactions <- NULL
  build$fastcc_epsilon <- NULL
  build$parent_prepruning <- "none"
  build$parent_algorithm_contract <- paste(
    "original CORDA dependency assessment on the complete constrained",
    "parent GEM; no FASTCC reaction deletion"
  )
  build$input_fastcore_epsilon_ignored_for_corda <-
    as.numeric(args$fastcore_epsilon %||% NA_real_)
  model$build_params <- build
  model$corda_parent_contract <- list(
    prepruning = "none",
    medium_constraints_applied = TRUE,
    forbidden_roles_blocked = c("demand", "sink", "artificial_support"),
    fastcore_epsilon_used = FALSE
  )
  model
}

.rc_complete_celltype_medium_corda_like_gem <-
  .rc_complete_celltype_medium_corda_gem
