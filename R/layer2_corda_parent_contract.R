# Complete medium-constrained parent model for original MATLAB CORDA2.

.rc_corda_parent <- function(
    gem, medium_table = NULL, condition = NULL,
    forbidden_roles = c("demand", "sink", "artificial_support"),
    solver = "highs", time_limit = 300) {
  parent <- rc_build_full_gem(
    gem, medium_table = medium_table, condition = condition
  )
  validated <- rc_validate_gem(parent)
  parent$corda_parent_original_lb <- validated$lb
  parent$corda_parent_original_ub <- validated$ub
  parent$corda_parent_n_reactions <- length(validated$reactions)
  parent$corda_parent_n_metabolites <- nrow(validated$S)
  parent$corda_parent_n_open_reactions <- sum(
    validated$lb < 0 | validated$ub > 0
  )
  parent$corda_parent_prepruning <- "none"
  parent$corda_parent_role_blocking <- "none"
  parent$corda_parent_feasibility_precheck <- FALSE
  parent$corda_parent_contract <- paste(
    "apply the requested medium to the complete GEM, then pass that model",
    "directly to original CORDA2 without FASTCC or reaction-role pruning"
  )
  parent$corda_ignored_fastcore_forbidden_roles <-
    unique(as.character(forbidden_roles))
  parent$corda_parent_solver <- as.character(solver)
  parent$corda_parent_time_limit <- as.numeric(time_limit)
  parent
}

.rc_corda_attach_parent_contract <- function(
    model, parent, fastcore_epsilon = NA_real_,
    forbidden_roles = c("demand", "sink", "artificial_support")) {
  if (!is.list(model$build_params)) {
    stop("CORDA2 build parameters are unavailable.", call. = FALSE)
  }
  build <- model$build_params
  build$n_parent_reactions <- as.integer(
    parent$corda_parent_n_reactions %||% NA_integer_
  )
  build$n_parent_metabolites <- as.integer(
    parent$corda_parent_n_metabolites %||% NA_integer_
  )
  build$n_parent_open_reactions <- as.integer(
    parent$corda_parent_n_open_reactions %||% NA_integer_
  )
  build$n_fastcc_consistent_parent_reactions <- NULL
  build$n_fastcc_inconsistent_parent_reactions <- NULL
  build$fastcc_epsilon <- NULL
  build$parent_prepruning <- "none"
  build$parent_role_blocking <- "none"
  build$parent_feasibility_precheck <- FALSE
  build$parent_algorithm_contract <- paste(
    "original MATLAB CORDA2 on the complete medium-constrained input GEM"
  )
  build$input_fastcore_epsilon_ignored_for_corda2 <-
    as.numeric(fastcore_epsilon)
  build$input_fastcore_forbidden_roles_ignored_for_corda2 <-
    unique(as.character(forbidden_roles))
  model$build_params <- build
  model$corda_parent_contract <- list(
    algorithm = "original_matlab_corda2",
    prepruning = "none",
    role_blocking = "none",
    feasibility_precheck = FALSE,
    n_reactions = build$n_parent_reactions,
    n_metabolites = build$n_parent_metabolites,
    n_open_reactions = build$n_parent_open_reactions,
    medium_constraints_applied = TRUE,
    fastcore_epsilon_used = FALSE,
    fastcore_forbidden_roles_used = FALSE
  )
  model$corda_noise_contract <- NULL
  model
}
