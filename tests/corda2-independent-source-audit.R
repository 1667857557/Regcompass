suppressPackageStartupMessages({
  library(Matrix)
  library(highs)
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
.rc_corda_empty_task_table <- function() data.frame()
.rc_lp_status <- function(message = "", code = NA_integer_) {
  text <- tolower(paste(message, collapse = " "))
  if (grepl("infeasible", text)) return("infeasible")
  if (grepl("unbounded", text)) return("unbounded")
  if (grepl("time|limit", text)) return("time_limit")
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
  if (anyNA(lb) || anyNA(ub) || any(lb > ub)) stop("invalid GEM bounds")
  list(S = S, lb = lb, ub = ub, reactions = reactions)
}
rc_solve_lp <- function(obj, A, lhs, rhs, lb, ub,
                        solver = "highs", time_limit = Inf) {
  answer <- highs::highs_solve(
    L = as.numeric(obj), lower = as.numeric(lb), upper = as.numeric(ub),
    A = A, lhs = as.numeric(lhs), rhs = as.numeric(rhs), maximum = FALSE,
    control = highs::highs_control(
      threads = 1L,
      time_limit = as.numeric(time_limit),
      log_to_console = FALSE,
      output_flag = FALSE,
      solver = "simplex",
      primal_feasibility_tolerance = 1e-7
    )
  )
  list(
    status = .rc_lp_status(answer$status_message, answer$status),
    solution = as.numeric(answer$primal_solution),
    objective = as.numeric(answer$objective_value),
    solver_message = answer$status_message
  )
}
rc_parallel_config <- function(...) list(workers = 1L)
rc_parallel_lapply <- function(X, FUN, BPPARAM = NULL, ...) lapply(X, FUN, ...)
.rc_layer2_task_bpparam <- function() FALSE

source("R/layer2_corda_evidence.R")
source("R/layer2_corda_lp.R")
source("R/layer2_corda_direction_contract.R")
source("R/layer2_corda2_algorithm.R")
source("R/layer2_corda2_algorithm_build.R")
source("R/layer2_corda2_options_contract.R")
source("R/layer2_corda_runtime.R")

oracle <- utils::read.delim(
  "tests/corda2_independent_oracle.tsv",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  colClasses = "character"
)
value <- function(case, key) {
  hit <- oracle$value[oracle$case == case & oracle$key == key]
  stopifnot(length(hit) == 1L)
  hit[[1L]]
}

make_gem <- function(metabolites, reactions, entries,
                     lb = NULL, ub = NULL) {
  S <- Matrix::Matrix(
    0,
    nrow = length(metabolites),
    ncol = length(reactions),
    sparse = TRUE,
    dimnames = list(metabolites, reactions)
  )
  for (entry in entries) S[entry[[1L]], entry[[2L]]] <- entry[[3L]]
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

# Constructor signature, exact source names and defaults.
defaults <- .rc_layer2_corda_options(list(model_completion = "corda2"))
stopifnot(
  identical(
    defaults$python_constructor_signature,
    strsplit(value("signature", "parameter_order"), ";", fixed = TRUE)[[1L]]
  ),
  is.null(defaults$met_prod),
  identical(defaults$n, as.integer(value("signature", "n"))),
  identical(defaults$penalty_factor,
            as.numeric(value("signature", "penalty_factor"))),
  identical(defaults$support,
            as.integer(value("signature", "support"))),
  identical(value("signature", "has_time_limit"), "false")
)
custom <- .rc_layer2_corda_options(list(
  model_completion = "corda2",
  corda2_args = list(
    met_prod = NULL, n = 1L, penalty_factor = 0.5, support = 2L
  )
))
stopifnot(
  identical(custom$n, 1L),
  identical(custom$penalty_factor, 0.5),
  identical(custom$support, 2L)
)
stopifnot(inherits(try(.rc_layer2_corda_options(list(
  model_completion = "corda2",
  corda2_args = list(n = 1L),
  corda2_redundancies = 2L
)), silent = TRUE), "try-error"))
stopifnot(inherits(try(.rc_layer2_corda_options(list(
  model_completion = "corda2",
  corda2_args = list(met_prod = "atp_c")
)), silent = TRUE), "try-error"))

# Ordered COBRA bound-setter behavior, including transient errors.
bound_cases <- list(
  bounds_reversible = c(-1000, 1000),
  bounds_positive = c(2, 1000),
  bounds_negative = c(-1000, -2),
  bounds_tiny = c(-5e-8, 5e-8),
  bounds_transient_lower_error = c(-2e6, -1.5e6),
  bounds_transient_upper_error = c(1.5e6, 2e6)
)
for (case in names(bound_cases)) {
  bounds <- bound_cases[[case]]
  gem <- make_gem(
    "A", "R", list(),
    lb = c(R = bounds[[1L]]), ub = c(R = bounds[[2L]])
  )
  observed <- try(.rc_corda_split_model(
    gem,
    tolerance = .rc_corda2_solver_feasibility_tolerance("highs"),
    upper_bound = 1e6
  ), silent = TRUE)
  if (identical(value(case, "status"), "error")) {
    stopifnot(inherits(observed, "try-error"))
  } else {
    stopifnot(!inherits(observed, "try-error"))
    checks <- c(
      normalized_reaction_lb = "reaction_lb",
      normalized_reaction_ub = "reaction_ub"
    )
    for (field in names(checks)) {
      stopifnot(isTRUE(all.equal(
        observed[[field]][["R"]],
        as.numeric(value(case, checks[[field]])),
        tolerance = 0
      )))
    }
    stopifnot(
      isTRUE(all.equal(observed$lb[["R::forward"]],
                       as.numeric(value(case, "forward_lb")), tolerance = 0)),
      isTRUE(all.equal(observed$ub[["R::forward"]],
                       as.numeric(value(case, "forward_ub")), tolerance = 0)),
      isTRUE(all.equal(observed$lb[["R::reverse"]],
                       as.numeric(value(case, "reverse_lb")), tolerance = 0)),
      isTRUE(all.equal(observed$ub[["R::reverse"]],
                       as.numeric(value(case, "reverse_ub")), tolerance = 0))
    )
  }
}

# Penalty-factor routing on a nondegenerate network.
route_gem <- make_gem(
  c("A", "B", "C", "D"),
  c("SRC", "N", "M1", "M2", "M3", "SINK"),
  list(
    list("A", "SRC", 1),
    list("A", "N", -1), list("D", "N", 1),
    list("A", "M1", -1), list("B", "M1", 1),
    list("B", "M2", -1), list("C", "M2", 1),
    list("C", "M3", -1), list("D", "M3", 1),
    list("D", "SINK", -1)
  )
)
route_classes <- classes_from_confidence(c(
  SRC = 0L, N = -1L, M1 = 1L, M2 = 1L, M3 = 1L, SINK = 3L
))
for (penalty_factor in c(0.5, 100)) {
  options <- .rc_layer2_corda_options(list(
    model_completion = "corda2",
    corda2_args = list(n = 1L, penalty_factor = penalty_factor, support = 5L)
  ))
  options$feasibility_tolerance <-
    .rc_corda2_solver_feasibility_tolerance("highs")
  split <- .rc_corda_split_model(
    route_gem,
    tolerance = options$feasibility_tolerance,
    upper_bound = options$upper_bound
  )
  confidence <- .rc_corda2_directional_confidence(split, route_classes)
  engine <- .rc_corda_new_lp_engine(split, "highs", Inf)
  associated <- .rc_corda2_associated(
    engine, split, "SINK::forward", confidence, options,
    penalize_medium = TRUE, redundancies = TRUE,
    stage = "independent_penalty_factor_audit"
  )
  engine <- .rc_corda_release_lp_engine(associated$engine)
  key <- paste0("penalty_factor_", format(penalty_factor, trim = TRUE))
  expected <- value(key, "needed")
  expected <- if (nzchar(expected)) {
    strsplit(expected, ";", fixed = TRUE)[[1L]]
  } else character()
  stopifnot(
    identical(sort(unique(associated$needed)), sort(expected)),
    identical(
      unname(associated$redundancies[["SINK::forward"]]),
      as.integer(value(key, "redundancy"))
    )
  )
}

# n=0 and verified persistent HiGHS configuration.
n0_gem <- make_gem(
  "A", c("SRC", "SINK"),
  list(list("A", "SRC", 1), list("A", "SINK", -1))
)
n0_options <- .rc_layer2_corda_options(list(
  model_completion = "corda2",
  corda2_args = list(n = 0L, penalty_factor = 100, support = 5L)
))
n0_options$feasibility_tolerance <-
  .rc_corda2_solver_feasibility_tolerance("highs")
n0_split <- .rc_corda_split_model(
  n0_gem, tolerance = n0_options$feasibility_tolerance
)
n0_confidence <- .rc_corda2_directional_confidence(
  n0_split, classes_from_confidence(c(SRC = 1L, SINK = 3L))
)
engine <- .rc_corda_new_lp_engine(n0_split, "highs", Inf)
stopifnot(
  identical(engine$type, "highs_persistent_cpp"),
  isTRUE(engine$solver_configuration_verified),
  identical(as.integer(engine$verified_solver_options$threads), 1L),
  identical(as.character(engine$verified_solver_options$solver), "simplex"),
  isTRUE(all.equal(
    as.numeric(engine$verified_solver_options$primal_feasibility_tolerance),
    n0_split$tolerance,
    tolerance = 0
  )),
  is.infinite(engine$time_limit)
)
n0 <- .rc_corda2_associated(
  engine, n0_split, "SINK::forward", n0_confidence, n0_options,
  penalize_medium = TRUE, redundancies = TRUE,
  stage = "independent_n_zero_audit"
)
engine <- .rc_corda_release_lp_engine(n0$engine)
stopifnot(
  !length(n0$needed),
  identical(unname(n0$redundancies[["SINK::forward"]]), 0L),
  identical(value("n_zero", "needed"), ""),
  identical(as.integer(value("n_zero", "redundancy")), 0L)
)

# The original model builder, not an override, passes Inf to CORDA2 build().
model_source <- paste(readLines("R/layer2_corda_model.R", warn = FALSE),
                      collapse = "\n")
stopifnot(grepl(
  "reconstruction <- .rc_corda_build_three_stage\\([\\s\\S]*time_limit = Inf",
  model_source,
  perl = TRUE
))

cat("Independent R/Python CORDA2 source audit passed\n")
