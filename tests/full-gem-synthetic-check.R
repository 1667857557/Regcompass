suppressPackageStartupMessages({
  library(Matrix)
  library(highs)
})

`%||%` <- function(x, y) if (is.null(x)) y else x
.rc_as_dgCMatrix <- function(x) methods::as(x, "dgCMatrix")
.rc_layer2_completion_context <- new.env(parent = emptyenv())

.rc_lp_status <- function(message = "", code = NA_integer_) {
  text <- tolower(paste(message, collapse = " "))
  if (grepl("infeasible", text)) return("infeasible")
  if (grepl("unbounded", text)) return("unbounded")
  if (grepl("time|limit", text)) return("time_limit")
  if (grepl("optimal", text)) return("optimal")
  if (is.finite(code) && as.integer(code) == 0L) return("optimal")
  "error"
}

rc_solve_lp <- function(obj, A, lhs, rhs, lb, ub,
                        solver = "highs", time_limit = Inf) {
  stopifnot(identical(solver, "highs"))
  answer <- highs::highs_solve(
    L = as.numeric(obj), lower = as.numeric(lb), upper = as.numeric(ub),
    A = A, lhs = as.numeric(lhs), rhs = as.numeric(rhs), maximum = FALSE,
    control = highs::highs_control(
      log_to_console = FALSE, output_flag = FALSE,
      threads = 1L, solver = "simplex",
      primal_feasibility_tolerance = 1e-7,
      time_limit = as.numeric(time_limit)
    )
  )
  list(
    status = .rc_lp_status(answer$status_message, answer$status),
    solution = as.numeric(answer$primal_solution),
    objective = as.numeric(answer$objective_value),
    solver_message = answer$status_message
  )
}

rc_align_bound <- function(x, rxns, default, name, allow_partial = FALSE) {
  if (is.null(x)) {
    return(stats::setNames(rep(default, length(rxns)), rxns))
  }
  x <- as.numeric(x[rxns])
  if (length(x) != length(rxns) || anyNA(x)) {
    stop("invalid ", name, " bounds", call. = FALSE)
  }
  stats::setNames(x, rxns)
}

rc_validate_gem <- function(gem) {
  S <- .rc_as_dgCMatrix(gem$S)
  reactions <- colnames(S)
  metabolites <- rownames(S)
  stopifnot(
    length(reactions) == ncol(S), !anyNA(reactions), !anyDuplicated(reactions),
    length(metabolites) == nrow(S), !anyNA(metabolites),
    !anyDuplicated(metabolites)
  )
  lb <- as.numeric(gem$lb[reactions])
  ub <- as.numeric(gem$ub[reactions])
  names(lb) <- reactions
  names(ub) <- reactions
  if (anyNA(lb) || anyNA(ub) || any(lb > ub)) stop("invalid GEM bounds")
  list(
    S = S, lb = lb, ub = ub,
    reactions = reactions, metabolites = metabolites
  )
}

rc_annotate_reaction_roles <- function(gem, medium_table = NULL) gem

rc_apply_medium_constraints <- function(
    gem, medium_table, condition = NULL, strict = FALSE, ...) {
  validated <- rc_validate_gem(gem)
  selected <- medium_table
  if ("condition" %in% colnames(selected)) {
    keep <- selected$condition == "all"
    if (!is.null(condition)) keep <- keep | selected$condition == condition
    selected <- selected[keep, , drop = FALSE]
  }
  for (i in seq_len(nrow(selected))) {
    reaction <- as.character(selected$exchange_reaction_id[[i]])
    if (!reaction %in% validated$reactions) next
    gem$lb[[reaction]] <- max(gem$lb[[reaction]], selected$lb[[i]])
    gem$ub[[reaction]] <- min(gem$ub[[reaction]], selected$ub[[i]])
  }
  list(
    gem = gem,
    medium_diagnostics = data.frame(
      reaction_id = selected$exchange_reaction_id,
      stringsAsFactors = FALSE
    )
  )
}

.rc_normalize_medium_scenarios <- function(x) x

rc_compass_vmax_directional <- function(
    S, lb, ub, target_reaction,
    direction = c("forward", "reverse"),
    solver = "highs", time_limit = 60,
    flux_threshold = 1e-8) {
  direction <- match.arg(direction)
  reactions <- colnames(S)
  index <- match(target_reaction, reactions)
  stopifnot(!is.na(index))
  objective <- rep(0, length(reactions))
  objective[[index]] <- if (identical(direction, "forward")) -1 else 1
  answer <- rc_solve_lp(
    obj = objective,
    A = S,
    lhs = rep(0, nrow(S)),
    rhs = rep(0, nrow(S)),
    lb = lb[reactions],
    ub = ub[reactions],
    solver = solver,
    time_limit = time_limit
  )
  vmax <- if (identical(answer$status, "optimal")) {
    if (identical(direction, "forward")) {
      answer$solution[[index]]
    } else {
      -answer$solution[[index]]
    }
  } else {
    0
  }
  list(
    feasible = identical(answer$status, "optimal") &&
      is.finite(vmax) && vmax >= flux_threshold,
    vmax = max(0, as.numeric(vmax)),
    status = answer$status
  )
}

source("R/fastcore.R")
# full_gem.R intentionally keeps the public/cache-facing functions canonical,
# while COMPASS exchange-bound helpers live in this helper source. Standalone
# source-level regressions must therefore load the helper explicitly, matching
# the package Collate contract rather than duplicating helpers into full_gem.R.
source("R/compass_medium_semantics.R")
source("R/full_gem.R")
source("R/layer2_corda_parent_contract.R")
source("R/microcompass_vmax_cache.R")

S <- Matrix::sparseMatrix(
  i = c(1, 1, 2, 2, 3, 3, 4, 4),
  j = c(1, 2, 2, 3, 4, 5, 5, 6),
  x = c(1, -1, 1, -1, 1, -1, 1, -1),
  dims = c(4, 6),
  dimnames = list(
    c("A_c", "B_c", "C_c", "D_c"),
    c("EX_A", "R1", "EX_B", "EX_C", "R2", "EX_D")
  )
)
reactions <- colnames(S)
gem <- list(
  S = S,
  lb = stats::setNames(rep(0, length(reactions)), reactions),
  ub = stats::setNames(rep(1000, length(reactions)), reactions),
  reaction_meta = data.frame(
    reaction_id = reactions,
    role = c("exchange", "internal", "exchange",
             "exchange", "internal", "exchange"),
    role_source = "synthetic",
    stringsAsFactors = FALSE
  ),
  metabolite_meta = data.frame(
    metabolite_id = rownames(S), stringsAsFactors = FALSE
  ),
  model_info = list(species = "synthetic", version = "1")
)
medium <- data.frame(
  medium_scenario_id = "close_C_source",
  exchange_reaction_id = c("EX_A", "EX_B", "EX_C", "EX_D"),
  lb = c(0, 0, 0, 0),
  ub = c(1000, 1000, 0, 1000),
  available = c(TRUE, TRUE, FALSE, TRUE),
  condition = "all",
  stringsAsFactors = FALSE
)

full <- rc_build_full_gem(gem, medium_table = medium)
stopifnot(
  identical(colnames(full$S), reactions),
  identical(rownames(full$S), rownames(S)),
  identical(full$S, .rc_as_dgCMatrix(S)),
  identical(unname(full$ub[["EX_C"]]), 0),
  identical(
    full$build_params$strategy,
    "compass_medium_constrained_full_gem"
  ),
  identical(full$build_params$n_input_reactions, 6L),
  identical(full$build_params$n_reactions, 6L),
  identical(full$build_params$n_medium_removed_reactions, 0L),
  identical(full$build_params$medium_direct_reaction_deletion, FALSE),
  identical(
    full$build_params$medium_handling,
    "exchange_bounds_only_no_reaction_deletion"
  )
)

r1_vmax <- rc_compass_vmax_directional(
  full$S, full$lb, full$ub, "R1", "forward",
  solver = "highs", flux_threshold = 1e-8
)
r2_vmax <- rc_compass_vmax_directional(
  full$S, full$lb, full$ub, "R2", "forward",
  solver = "highs", flux_threshold = 1e-8
)
r2_step2 <- .rc_compass_step2_from_vmax_directional(
  S = full$S,
  lb = full$lb,
  ub = full$ub,
  target_reaction = "R2",
  penalties = stats::setNames(rep(0, length(reactions)), reactions),
  vmax_result = r2_vmax,
  target_direction = "forward",
  solver = "highs",
  flux_threshold = 1e-8
)
stopifnot(
  isTRUE(r1_vmax$feasible),
  r1_vmax$vmax > 0,
  !isTRUE(r2_vmax$feasible),
  r2_vmax$vmax < 1e-8,
  identical(r2_step2$feasible, FALSE),
  identical(r2_step2$step2_status, "not_run")
)

dirs <- data.frame(
  reaction_id = c("R1", "R2"),
  target_direction = "forward",
  stringsAsFactors = FALSE
)
cache_dir <- tempfile("full-gem-ci-")
cache <- rc_build_full_gem_cache(
  gem = gem,
  dirs = dirs,
  medium_scenarios = medium,
  cache_dir = cache_dir,
  solver = "highs",
  time_limit = 60,
  flux_consistency_epsilon = 1e-8
)
summary <- attr(cache, "summary")
stopifnot(
  length(cache) == 2L,
  setequal(vapply(cache, `[[`, character(1), "reaction_id"), c("R1", "R2")),
  identical(summary$n_input_reactions, 6L),
  identical(summary$n_reactions, 6L),
  identical(summary$n_medium_removed_reactions, 0L),
  identical(summary$medium_applied, TRUE),
  identical(
    summary$medium_handling,
    "exchange_bounds_only_no_reaction_deletion"
  ),
  is.character(summary$medium_fingerprint),
  nzchar(summary$medium_fingerprint),
  identical(attr(cache, "completion_method"), "none"),
  identical(attr(cache, "fastcore_executed"), FALSE),
  identical(attr(cache, "corda2_executed"), FALSE),
  identical(
    attr(cache, "medium_handling"),
    "exchange_bounds_only_no_reaction_deletion"
  )
)

medium_open <- medium
medium_open$ub[medium_open$exchange_reaction_id == "EX_C"] <- 1000
medium_open$available[medium_open$exchange_reaction_id == "EX_C"] <- TRUE
cache_open <- rc_build_full_gem_cache(
  gem = gem,
  dirs = dirs,
  medium_scenarios = medium_open,
  cache_dir = cache_dir,
  solver = "highs",
  time_limit = 60,
  flux_consistency_epsilon = 1e-8
)
summary_open <- attr(cache_open, "summary")
full_open <- readRDS(summary_open$file[[1L]])
r2_open <- rc_compass_vmax_directional(
  full_open$S, full_open$lb, full_open$ub,
  "R2", "forward", solver = "highs", flux_threshold = 1e-8
)
stopifnot(
  length(cache_open) == 2L,
  identical(summary_open$n_reactions, 6L),
  identical(summary_open$n_medium_removed_reactions, 0L),
  !identical(summary$file, summary_open$file),
  !identical(summary$medium_fingerprint, summary_open$medium_fingerprint),
  isTRUE(r2_open$feasible)
)

fastcore_parent <- .rc_fastcore_parent(
  gem,
  medium_table = medium,
  forbidden_roles = character(),
  solver = "highs",
  time_limit = 60,
  fastcore_epsilon = 1e-4
)
corda_parent <- .rc_corda_parent(
  gem,
  medium_table = medium,
  forbidden_roles = character(),
  solver = "highs",
  time_limit = 60
)
stopifnot(
  identical(colnames(fastcore_parent$S), reactions),
  identical(colnames(corda_parent$S), reactions),
  identical(corda_parent$corda_parent_prepruning, "none"),
  identical(corda_parent$corda_parent_role_blocking, "none")
)

for (invalid in list(
  list(model_completion = "fastcore"),
  list(model_completion = "corda2"),
  list(fastcore_epsilon = 1e-4),
  list(corda2_args = list(MCxNCthresh = 2))
)) {
  error <- try(.rc_validate_full_gem_model_params(invalid), silent = TRUE)
  stopifnot(inherits(error, "try-error"))
}
valid <- .rc_validate_full_gem_model_params(
  list(model_completion = "none", completion_time_limit = 60)
)
stopifnot(identical(valid$model_completion, "none"))

message("COMPASS medium-bounds synthetic contract passed.")
