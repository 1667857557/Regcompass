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
  if (length(lb) != length(reactions) || length(ub) != length(reactions)) {
    stop("synthetic GEM bound length mismatch")
  }
  names(lb) <- reactions
  names(ub) <- reactions
  if (anyNA(lb) || anyNA(ub) || any(lb > ub)) stop("invalid GEM bounds")
  list(S = S, lb = lb, ub = ub, reactions = reactions)
}
rc_solve_lp <- function(obj, A, lhs, rhs, lb, ub,
                        solver = "highs", time_limit = Inf) {
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
rc_parallel_config <- function(...) list(workers = 1L)

source("R/layer2_corda_meta_module_contract.R")
source("R/layer2_corda_evidence.R")
source("R/layer2_corda_lp.R")
source("R/layer2_corda_paper_contract.R")
source("R/layer2_corda_direction_contract.R")
source("R/layer2_corda_target_isolation.R")
source("R/layer2_corda2_algorithm.R")
source("R/layer2_corda2_algorithm_build.R")
source("R/layer2_corda2_options_contract.R")

make_gem <- function(metabolites, reactions, entries, lb = NULL, ub = NULL) {
  S <- Matrix::sparseMatrix(
    i = integer(), j = integer(), x = numeric(),
    dims = c(length(metabolites), length(reactions)),
    dimnames = list(metabolites, reactions),
    giveCsparse = TRUE
  )
  for (entry in entries) S[entry[[1L]], entry[[2L]]] <- entry[[3L]]
  if (is.null(lb)) lb <- stats::setNames(rep(0, length(reactions)), reactions)
  if (is.null(ub)) ub <- stats::setNames(rep(1000, length(reactions)), reactions)
  stopifnot(
    identical(colnames(S), reactions),
    length(lb) == length(reactions),
    length(ub) == length(reactions)
  )
  list(S = S, lb = lb, ub = ub)
}
classes <- function(hc, mc, nc, ot = character()) {
  reactions <- unique(c(hc, mc, nc, ot))
  label <- stats::setNames(rep("OT", length(reactions)), reactions)
  label[hc] <- "HC"
  label[mc] <- "MC_module"
  label[nc] <- "NC"
  list(
    hc = hc, mc_module = mc, mc_evidence = character(), mc = mc,
    nc = nc, ot = ot, confidence = label, initial_confidence = label
  )
}
inspect_gem <- function(label, gem) {
  cat(
    label, ": reactions=", paste(colnames(gem$S), collapse = ","),
    "; ncol=", ncol(gem$S),
    "; lb=", paste(names(gem$lb), collapse = ","),
    "; ub=", paste(names(gem$ub), collapse = ","), "\n",
    sep = ""
  )
}

options <- .rc_layer2_corda_options(list(model_completion = "corda2"))
stopifnot(
  identical(options$MCxNCthresh, 2),
  identical(options$constraint, 1),
  identical(options$constrainby, "val"),
  identical(options$om, 1e4),
  identical(options$ci, 0.01),
  identical(options$flux_threshold, 1e-7),
  identical(options$baseline_cost, 1e-3),
  identical(options$time_limit, Inf),
  identical(options$stage_order,
            c("step1_HC_dependencies", "step2_1_MC_NC_dependencies",
              "step2_2_MC_feasibility", "step3_HC_OT_dependencies"))
)

# 1. Reversible split keeps one forward and one reverse directional variable and
# closes the opposite copy under target assessment. Validate semantic direction
# from the direction table rather than hard-coding an internal variable label.
gem <- make_gem(
  c("A", "B"), c("REV", "IRR"),
  list(c(1, 1, -1), c(2, 1, 1), c(2, 2, -1)),
  lb = c(REV = -1000, IRR = 0),
  ub = c(REV = 1000, IRR = 1000)
)
inspect_gem("reversible", gem)
split <- .rc_corda2_split_original(gem)
forward_id <- split$direction_table$variable_id[
  split$direction_table$reaction_id == "REV" &
    split$direction_table$direction == "forward"
]
reverse_id <- split$direction_table$variable_id[
  split$direction_table$reaction_id == "REV" &
    split$direction_table$direction == "reverse"
]
stopifnot(
  length(forward_id) == 1L,
  length(reverse_id) == 1L,
  identical(.rc_corda2_opposite_variable(split, forward_id), reverse_id),
  identical(.rc_corda2_opposite_variable(split, reverse_id), forward_id)
)
constrained <- .rc_corda2_constrain_target(
  .rc_corda_new_lp_engine(split, "highs", Inf),
  split, forward_id, options
)
stopifnot(
  constrained$upper[[reverse_id]] == 0,
  constrained$lower[[forward_id]] >= 0,
  identical(constrained$answer$status, "optimal")
)
.rc_corda_release_lp_engine(constrained$engine)

# 2. Step 1 HC dependency promotion and Step 2.1/2.2/3 state transitions are
# exercised by a compact source-to-sink network.
gem <- make_gem(
  c("A", "B", "C"),
  c("UP", "HC", "MC", "NC", "OUT"),
  list(
    c(1, 1, 1),
    c(1, 2, -1), c(2, 2, 1),
    c(2, 3, -1), c(3, 3, 1),
    c(2, 4, -1), c(3, 4, 1),
    c(3, 5, -1)
  ),
  lb = stats::setNames(rep(0, 5), c("UP", "HC", "MC", "NC", "OUT")),
  ub = stats::setNames(rep(1000, 5), c("UP", "HC", "MC", "NC", "OUT"))
)
reconstruction <- .rc_corda_build_three_stage(
  split = .rc_corda2_split_original(gem),
  classes = classes(hc = c("UP", "HC", "OUT"), mc = "MC", nc = "NC"),
  options = options,
  solver = "highs",
  time_limit = Inf
)
stopifnot(
  is.list(reconstruction),
  identical(reconstruction$stage_order, options$stage_order),
  is.matrix(reconstruction$HCtoMC),
  is.matrix(reconstruction$HCtoNC),
  is.matrix(reconstruction$MCtoNC),
  all(reconstruction$final_confidence %in% 0:3),
  reconstruction$solver_performance$n_solves >= 1L
)

# 3. The canonical option object contains all original user-facing controls.
for (name in c("MCxNCthresh", "constraint", "constrainby", "om", "ci")) {
  stopifnot(name %in% names(options))
}
cat("Original CORDA2 synthetic source checks passed.\n")
