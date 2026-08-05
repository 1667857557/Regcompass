suppressPackageStartupMessages({
  library(Matrix)
  library(highs)
})

`%||%` <- function(x, y) if (is.null(x)) y else x
.rc_as_dgCMatrix <- function(x) methods::as(x, "dgCMatrix")
.rc_layer2_completion_context <- new.env(parent = emptyenv())
.rc_layer2_completion_context$solver <- "highs"
.rc_layer2_completion_context$completion_time_limit <- 60
.rc_layer2_completion_context$flux_threshold <- 1e-8

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

source("R/fastcore.R")
source("R/full_gem.R")

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

parent <- rc_build_full_gem(gem)
stopifnot(
  identical(colnames(parent$S), reactions),
  identical(parent$build_params$strategy, "full_gem"),
  identical(parent$target_status, "not_prechecked"),
  is.na(parent$build_params$n_flux_inconsistent_reactions)
)

pruned <- rc_build_full_gem(
  gem,
  medium_table = medium,
  prune_flux_inconsistent = TRUE,
  solver = "highs",
  time_limit = 60,
  flux_consistency_epsilon = 1e-8
)
stopifnot(
  setequal(colnames(pruned$S), c("EX_A", "R1", "EX_B")),
  !any(c("EX_C", "R2", "EX_D") %in% colnames(pruned$S)),
  identical(
    pruned$build_params$strategy,
    "medium_flux_consistency_pruned_full_gem"
  ),
  identical(pruned$build_params$n_input_reactions, 6L),
  identical(pruned$build_params$n_flux_inconsistent_reactions, 3L),
  identical(pruned$build_params$fastcore_executed, FALSE),
  identical(pruned$build_params$corda2_executed, FALSE),
  identical(pruned$build_params$reaction_evidence_used, FALSE)
)

dirs <- data.frame(
  reaction_id = c("R1", "R2"),
  target_direction = "forward",
  stringsAsFactors = FALSE
)
cache <- rc_build_full_gem_cache(
  gem = gem,
  dirs = dirs,
  medium_scenarios = medium,
  cache_dir = tempfile("full-gem-ci-"),
  solver = "highs",
  time_limit = 60,
  flux_consistency_epsilon = 1e-8
)
summary <- attr(cache, "summary")
stopifnot(
  length(cache) == 1L,
  identical(cache[[1L]]$reaction_id, "R1"),
  identical(summary$n_input_reactions, 6L),
  identical(summary$n_reactions, 3L),
  identical(summary$n_flux_inconsistent_reactions, 3L),
  identical(summary$medium_applied, TRUE),
  identical(summary$solver, "highs"),
  identical(summary$completion_time_limit, 60),
  identical(attr(cache, "completion_method"), "none"),
  identical(attr(cache, "fastcore_executed"), FALSE),
  identical(attr(cache, "corda2_executed"), FALSE)
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
stopifnot(invisible(.rc_validate_full_gem_model_params(
  list(model_completion = "none", completion_time_limit = 60)
)))

message("Full-GEM medium-pruning synthetic contract passed.")
