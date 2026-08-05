suppressPackageStartupMessages({
  library(Matrix)
  library(highs)
  library(BiocParallel)
})

`%||%` <- function(x, y) if (is.null(x)) y else x
.rc_as_dgCMatrix <- function(x) methods::as(x, "dgCMatrix")
.rc_bind_frames_fill <- function(values) {
  values <- values[vapply(values, is.data.frame, logical(1))]
  values <- values[vapply(values, nrow, integer(1)) > 0L]
  if (!length(values)) return(data.frame())
  columns <- unique(unlist(lapply(values, colnames), use.names = FALSE))
  values <- lapply(values, function(value) {
    missing <- setdiff(columns, colnames(value))
    for (name in missing) value[[name]] <- NA
    value[, columns, drop = FALSE]
  })
  answer <- do.call(rbind, values)
  rownames(answer) <- NULL
  answer
}
.rc_lp_status <- function(message = "", code = NA_integer_) {
  text <- tolower(paste(message, collapse = " "))
  if (grepl("infeasible", text)) return("infeasible")
  if (grepl("unbounded", text)) return("unbounded")
  if (grepl("optimal", text)) return("optimal")
  if (is.finite(code) && as.integer(code) == 0L) return("optimal")
  "error"
}
rc_validate_gem <- function(gem) {
  S <- .rc_as_dgCMatrix(gem$S)
  reactions <- colnames(S)
  lb <- as.numeric(gem$lb[reactions])
  ub <- as.numeric(gem$ub[reactions])
  names(lb) <- names(ub) <- reactions
  list(S = S, lb = lb, ub = ub, reactions = reactions)
}
rc_solve_lp <- function(obj, A, lhs, rhs, lb, ub,
                        solver = "highs", time_limit = Inf) {
  control <- list(log_to_console = FALSE, threads = 1L, solver = "simplex")
  if (is.finite(time_limit)) control$time_limit <- time_limit
  answer <- highs::highs_solve(
    L = as.numeric(obj), lower = as.numeric(lb), upper = as.numeric(ub),
    A = A, lhs = as.numeric(lhs), rhs = as.numeric(rhs), maximum = FALSE,
    control = do.call(highs::highs_control, control)
  )
  list(
    status = .rc_lp_status(answer$status_message, answer$status),
    solution = as.numeric(answer$primal_solution),
    objective = as.numeric(answer$objective_value),
    solver_message = answer$status_message
  )
}
rc_parallel_config <- function(workers = NULL, backend = "auto") {
  list(workers = 1L)
}
rc_parallel_lapply <- function(X, FUN, BPPARAM = NULL, ...) {
  lapply(X, FUN, ...)
}
.rc_layer2_task_bpparam <- function() FALSE
.rc_subset_gem <- function(gem, reactions) gem
rc_prepare_directional_targets <- function(...) data.frame()
.rc_directional_feasibility <- function(...) data.frame()
rc_build_full_gem <- function(gem, ...) gem
.rc_fastcore_parent <- function(...) stop("not used")
rc_export_microcompass <- function(...) invisible(TRUE)

source("R/layer2_corda_evidence.R")
source("R/layer2_corda_lp.R")
source("R/layer2_corda_paper_contract.R")
source("R/layer2_corda_direction_contract.R")
source("R/layer2_corda_model.R")
source("R/layer2_corda_output_contract.R")
source("R/layer2_corda_target_contract.R")
source("R/layer2_corda_parent_contract.R")
source("R/layer2_corda2_algorithm.R")
source("R/layer2_corda2_algorithm_build.R")
source("R/layer2_corda2_algorithm_integration.R")
source("R/layer2_corda2_options_contract.R")
source("R/layer2_corda2_correction_contract.R")
source("R/layer2_corda2_output_contract.R")

oracle <- utils::read.delim(
  "tests/corda2_python_oracle.tsv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
value <- function(case, key) {
  hit <- oracle$value[oracle$case == case & oracle$key == key]
  stopifnot(length(hit) == 1L)
  hit[[1L]]
}

make_gem <- function(metabolites, reactions, entries, lb = NULL, ub = NULL) {
  S <- Matrix::Matrix(
    0, nrow = length(metabolites), ncol = length(reactions), sparse = TRUE,
    dimnames = list(metabolites, reactions)
  )
  for (entry in entries) {
    S[entry[[1L]], entry[[2L]]] <- entry[[3L]]
  }
  if (is.null(lb)) lb <- stats::setNames(rep(0, length(reactions)), reactions)
  if (is.null(ub)) ub <- stats::setNames(rep(1000, length(reactions)), reactions)
  list(S = S, lb = lb, ub = ub)
}

classes_from_confidence <- function(confidence) {
  reactions <- names(confidence)
  label <- stats::setNames(rep("OT", length(reactions)), reactions)
  label[confidence == -1L] <- "NC"
  label[confidence == 1L] <- "MC_evidence"
  label[confidence == 2L] <- "MC_module"
  label[confidence == 3L] <- "HC"
  list(
    hc = reactions[confidence == 3L],
    mc_module = reactions[confidence == 2L],
    mc_evidence = reactions[confidence == 1L],
    mc = reactions[confidence %in% c(1L, 2L)],
    nc = reactions[confidence == -1L],
    ot = reactions[confidence == 0L],
    confidence = label,
    initial_confidence = label
  )
}

options_for <- function(n = 3L, support = 5L, penalty = 100) {
  .rc_layer2_corda_options(list(
    model_completion = "corda2",
    corda2_redundancies = n,
    corda2_support = support,
    corda2_penalty_factor = penalty
  ))
}

run_build <- function(gem, confidence, n = 3L, support = 5L, penalty = 100) {
  options <- options_for(n, support, penalty)
  options$feasibility_tolerance <-
    .rc_corda2_solver_feasibility_tolerance("highs")
  split <- .rc_corda_split_model(
    gem,
    tolerance = options$feasibility_tolerance,
    upper_bound = options$upper_bound
  )
  result <- .rc_corda_build_three_stage(
    split, classes_from_confidence(confidence), options,
    solver = "highs", time_limit = 60
  )
  list(result = result, split = split, options = options)
}

compare_oracle_confidence <- function(case, result) {
  rows <- oracle[oracle$case == case & grepl("::", oracle$key, fixed = TRUE), ]
  observed <- result$final_directional_confidence[rows$key]
  stopifnot(
    !anyNA(observed),
    identical(as.integer(observed), as.integer(rows$value))
  )
  expected_included <- value(case, "included")
  expected_included <- if (nzchar(expected_included)) {
    strsplit(expected_included, ";", fixed = TRUE)[[1L]]
  } else {
    character()
  }
  stopifnot(setequal(result$included, expected_included))
}

# Source constants and active-solver feasibility tolerance.
stopifnot(
  identical(as.numeric(value("constants", "UPPER")), 1e6),
  identical(as.numeric(value("constants", "CI")), 1.01),
  identical(as.numeric(value("constants", "tflux")), 1),
  identical(as.integer(value("constants", "default_n")), 3L),
  identical(as.numeric(value("constants", "default_penalty_factor")), 100),
  identical(as.integer(value("constants", "default_support")), 5L),
  isTRUE(all.equal(
    as.numeric(value("constants", "feasibility_tolerance")),
    .rc_corda2_solver_feasibility_tolerance("highs"),
    tolerance = 0
  ))
)

# Upstream redundancy example: two alternative routes are recovered.
association_gem <- make_gem(
  c("A", "B", "C"),
  c("R1", "R2", "SRC_A", "SRC_B", "SINK_C"),
  list(
    list("A", "R1", -1), list("C", "R1", 1),
    list("B", "R2", -1), list("C", "R2", 1),
    list("A", "SRC_A", 1), list("B", "SRC_B", 1),
    list("C", "SINK_C", -1)
  )
)
association_conf <- c(R1 = 1L, R2 = 2L, SRC_A = 1L, SRC_B = 1L, SINK_C = 3L)
association_options <- options_for(n = 3L)
association_options$feasibility_tolerance <-
  .rc_corda2_solver_feasibility_tolerance("highs")
association_split <- .rc_corda_split_model(
  association_gem,
  tolerance = association_options$feasibility_tolerance
)
association_directional <- .rc_corda2_directional_confidence(
  association_split, classes_from_confidence(association_conf)
)
association_engine <- .rc_corda_new_lp_engine(
  association_split, "highs", 60
)
association <- .rc_corda2_associated(
  association_engine,
  association_split,
  "SINK_C::forward",
  association_directional,
  association_options,
  penalize_medium = TRUE,
  redundancies = TRUE,
  stage = "oracle_association"
)
expected_needed <- strsplit(
  value("association", "needed"), ";", fixed = TRUE
)[[1L]]
stopifnot(
  setequal(association$needed, expected_needed),
  identical(
    unname(association$redundancies[["SINK_C::forward"]]),
    as.integer(value("association", "redundancy"))
  )
)

# The opposite reversible variable is intentionally left open, exactly as the
# pinned Python source does.
self_cycle_gem <- make_gem(
  c("A", "B"), "REV",
  list(list("A", "REV", -1), list("B", "REV", 1)),
  lb = c(REV = -1000), ub = c(REV = 1000)
)
self_options <- options_for(n = 1L)
self_options$feasibility_tolerance <-
  .rc_corda2_solver_feasibility_tolerance("highs")
self_split <- .rc_corda_split_model(
  self_cycle_gem, tolerance = self_options$feasibility_tolerance
)
self_bounds <- .rc_corda_target_bounds(
  self_split, "REV::forward", epsilon = 1
)
stopifnot(
  identical(self_bounds$opposite_direction_blocked, character()),
  self_bounds$upper[["REV::reverse"]] == 1e6
)
self_answer <- rc_solve_lp(
  obj = c(0, 0), A = self_split$S,
  lhs = rep(0, nrow(self_split$S)), rhs = rep(0, nrow(self_split$S)),
  lb = self_bounds$lower, ub = self_bounds$upper,
  solver = "highs", time_limit = 60
)
stopifnot(
  identical(self_answer$status, value("self_cycle", "status")),
  self_answer$solution[[1L]] > self_split$tolerance,
  self_answer$solution[[2L]] > self_split$tolerance
)

# Unknown/free completion in build iteration 3.
stage3_gem <- make_gem(
  "A", c("SRC", "H"),
  list(list("A", "SRC", 1), list("A", "H", -1))
)
stage3 <- run_build(stage3_gem, c(SRC = 0L, H = 3L))$result
compare_oracle_confidence("stage3_unknown", stage3)
stopifnot("SRC" %in% stage3$stage3_associated_ot)

# Absent-support counting plus the source's positive-minimization step.
support_gem <- make_gem(
  "X", c("N", "M1", "M2"),
  list(
    list("X", "N", 1),
    list("X", "M1", -1),
    list("X", "M2", -1)
  ),
  lb = c(N = 0, M1 = 2, M2 = 0),
  ub = c(N = 1000, M1 = 1000, M2 = 1000)
)
support <- run_build(
  support_gem, c(N = -1L, M1 = 2L, M2 = 2L),
  n = 1L, support = 2L
)$result
compare_oracle_confidence("absent_support", support)
stopifnot(
  support$stage2_nc_support_count[["N::forward"]] == 2L,
  "N" %in% support$stage2_promoted_nc,
  "M1" %in% support$stage2_promoted_mc,
  !"M2" %in% support$stage2_promoted_mc
)

# Exact source behavior in iteration 2b: minimize +1*v and promote only when
# the minimum objective is greater than tflux=1.
forced_gem <- make_gem(
  "A", c("SRC", "M"),
  list(list("A", "SRC", 1), list("A", "M", -1)),
  lb = c(SRC = 0, M = 2), ub = c(SRC = 1000, M = 1000)
)
forced <- run_build(forced_gem, c(SRC = 0L, M = 2L), n = 1L)$result
compare_oracle_confidence("positive_min_forced", forced)
stopifnot("M" %in% forced$stage2_promoted_mc)

free_gem <- make_gem(
  "A", c("SRC", "M"),
  list(list("A", "SRC", 1), list("A", "M", -1))
)
free <- run_build(free_gem, c(SRC = 0L, M = 2L), n = 1L)$result
compare_oracle_confidence("positive_min_free", free)
stopifnot(!"M" %in% free$stage2_promoted_mc)

# Python assigns both direction costs from the forward variable confidence.
penalty_split <- .rc_corda_split_model(
  self_cycle_gem,
  tolerance = .rc_corda2_solver_feasibility_tolerance("highs")
)
penalty <- .rc_corda2_penalties(
  penalty_split,
  c("REV::forward" = 3L, "REV::reverse" = -1L),
  penalize_medium = TRUE,
  penalty_factor = 100
)
stopifnot(
  penalty[["REV::forward"]] == 0,
  penalty[["REV::reverse"]] == 0
)

# The one-shot fallback preserves the same exact control flow and result on the
# nondegenerate oracle network.
original_api <- .rc_corda_highs_api_available
.rc_corda_highs_api_available <- function() FALSE
stage3_one_shot <- run_build(stage3_gem, c(SRC = 0L, H = 3L))$result
.rc_corda_highs_api_available <- original_api
stopifnot(
  identical(
    stage3_one_shot$final_directional_confidence,
    stage3$final_directional_confidence
  ),
  setequal(stage3_one_shot$included, stage3$included)
)

stopifnot(
  identical(
    stage3$algorithm,
    "resendislab_python_CORDA2_c02e06d_exact_semantics"
  ),
  identical(stage3$stage_update_policy, "python_serial_mutation_order"),
  all(vapply(stage3$execution, function(x) {
    identical(x$workers, 1L) && identical(x$target_parallelism, FALSE)
  }, logical(1)))
)

cat("R implementation matches the pinned Python CORDA2 oracle\n")
