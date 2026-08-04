# Preserve the original CORDA parent-model semantics without pre-pruning.

.rc_fastcore_parent_before_corda_contract <- .rc_fastcore_parent

.rc_corda_parent <- function(
    gem, medium_table = NULL, condition = NULL,
    forbidden_roles = c("demand", "sink", "artificial_support"),
    solver = "highs", time_limit = 300) {
  parent <- rc_build_full_gem(
    gem, medium_table = medium_table, condition = condition
  )
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
  parent$corda_parent_original_lb <- validated$lb
  parent$corda_parent_original_ub <- validated$ub
  parent$corda_parent_n_reactions <- length(validated$reactions)
  parent$corda_parent_n_metabolites <- nrow(validated$S)
  parent$corda_parent_n_open_reactions <- sum(
    validated$lb < 0 | validated$ub > 0
  )
  parent$corda_parent_prepruning <- "none"
  parent$corda_parent_role_blocking <- "none"
  parent$corda_parent_contract <- paste(
    "the complete input GEM is retained after applying only the requested",
    "medium bounds; no FASTCC pruning or reaction-role blocking is applied"
  )
  parent$corda_ignored_fastcore_forbidden_roles <-
    unique(as.character(forbidden_roles))
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

.rc_corda_medium_id <- function(medium_table) {
  if (is.null(medium_table) || !is.data.frame(medium_table) ||
      !nrow(medium_table) ||
      !"medium_scenario_id" %in% colnames(medium_table)) {
    return("default")
  }
  value <- unique(as.character(medium_table$medium_scenario_id))
  value <- value[!is.na(value) & nzchar(value)]
  if (length(value) != 1L) {
    stop(
      "CORDA model construction requires exactly one medium scenario.",
      call. = FALSE
    )
  }
  value
}

.rc_corda_noise_namespace <- function(cell_type, medium_table) {
  cell_type <- as.character(cell_type)
  if (length(cell_type) != 1L || is.na(cell_type) || !nzchar(cell_type)) {
    stop("CORDA noise namespace requires one non-empty cell type.",
         call. = FALSE)
  }
  medium <- .rc_corda_medium_id(medium_table)
  paste0(
    "celltype=", utils::URLencode(cell_type, reserved = TRUE),
    "::medium=", utils::URLencode(medium, reserved = TRUE)
  )
}

.rc_complete_celltype_medium_corda_gem_parent_base <-
  .rc_complete_celltype_medium_corda_gem

.rc_complete_celltype_medium_corda_gem <- function(...) {
  args <- list(...)
  corda_options <- args$corda_options
  if (!is.list(corda_options)) {
    stop("CORDA options are unavailable during model construction.",
         call. = FALSE)
  }
  noise_namespace <- .rc_corda_noise_namespace(
    cell_type = args$cell_type,
    medium_table = args$medium_table
  )
  corda_options$noise_namespace <- noise_namespace
  args$corda_options <- corda_options

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
  build$n_parent_reactions <- as.integer(
    model$corda_parent_n_reactions %||% NA_integer_
  )
  build$n_parent_metabolites <- as.integer(
    model$corda_parent_n_metabolites %||% NA_integer_
  )
  build$n_parent_open_reactions <- as.integer(
    model$corda_parent_n_open_reactions %||% NA_integer_
  )
  build$n_fastcc_consistent_parent_reactions <- NULL
  build$n_fastcc_inconsistent_parent_reactions <- NULL
  build$fastcc_epsilon <- NULL
  build$parent_prepruning <- "none"
  build$parent_role_blocking <- "none"
  build$parent_algorithm_contract <- paste(
    "original CORDA dependency assessment on the complete medium-constrained",
    "input GEM; no FASTCC deletion and no generic role-based reaction block"
  )
  build$input_fastcore_epsilon_ignored_for_corda <-
    as.numeric(args$fastcore_epsilon %||% NA_real_)
  build$input_fastcore_forbidden_roles_ignored_for_corda <-
    model$corda_ignored_fastcore_forbidden_roles %||% character()
  build$corda_noise_namespace <- noise_namespace
  build$corda_noise_task_key <-
    "namespace_x_stage_x_signed_target_x_repeat"
  model$build_params <- build
  model$corda_parent_contract <- list(
    prepruning = "none",
    role_blocking = "none",
    n_reactions = build$n_parent_reactions,
    n_metabolites = build$n_parent_metabolites,
    n_open_reactions = build$n_parent_open_reactions,
    medium_constraints_applied = TRUE,
    fastcore_epsilon_used = FALSE,
    fastcore_forbidden_roles_used = FALSE
  )
  model$corda_noise_contract <- list(
    seed = corda_options$seed,
    namespace = noise_namespace,
    task_key = "namespace_x_stage_x_signed_target_x_repeat",
    distribution = paste0("Uniform(0, ", corda_options$kappa, ")"),
    scheduling_invariant = TRUE,
    matlab_rng_bitwise_identity = FALSE
  )
  model
}

.rc_complete_celltype_medium_corda_like_gem <-
  .rc_complete_celltype_medium_corda_gem
