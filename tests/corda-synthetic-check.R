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
source("R/layer2_corda2_algorithm.R")
source("R/layer2_corda2_algorithm_build.R")
source("R/layer2_corda_serial_core.R")
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
  identical(options$output_bound, 1000)
)
custom <- .rc_layer2_corda_options(list(
  model_completion = "corda2",
  corda2_args = list(
    MCxNCthresh = 3, constraint = 25,
    constrainby = "perc", om = 1e5, ci = 0.02
  )
))
stopifnot(
  custom$MCxNCthresh == 3,
  custom$constraint == 25,
  custom$constrainby == "perc",
  custom$om == 1e5,
  custom$ci == 0.02
)
stopifnot(inherits(try(
  .rc_layer2_corda_options(list(
    model_completion = "corda2", corda2_redundancies = 3L
  )), silent = TRUE
), "try-error"))

handoff_layer1 <- list(
  reaction_support = data.frame(
    reaction_id = c("R1", "R2"),
    cell_type = c("epithelial", "epithelial"),
    rna_reaction_support = c(0.8, 0.2),
    multiome_reaction_support = c(0.9, 0.3),
    stringsAsFactors = FALSE
  )
)
handoff_meta_modules <- list(
  workflow_params = list(celltype_col = "broad_type"),
  merged_modules = list(
    celltype_col = "broad_type",
    merged_reaction_membership = data.frame(
      broad_type = "epithelial",
      reaction_id = "R1",
      stringsAsFactors = FALSE
    )
  )
)
handoff_evidence <- .rc_corda_reaction_evidence(
  layer1 = handoff_layer1,
  meta_modules = handoff_meta_modules,
  regulatory_weight = 0.20
)
stopifnot(
  handoff_evidence$merged_meta_module_member[
    handoff_evidence$reaction_id == "R1"
  ],
  !handoff_evidence$merged_meta_module_member[
    handoff_evidence$reaction_id == "R2"
  ],
  identical(
    unique(handoff_evidence$evidence_source),
    "legacy_reaction_support_tables"
  )
)

reversible <- make_gem(
  c("A", "B"), c("REV", "IRR"),
  list(
    list("A", "REV", -1), list("B", "REV", 1),
    list("B", "IRR", -1)
  ),
  lb = c(REV = -1000, IRR = 0),
  ub = c(REV = 1000, IRR = 1000)
)
inspect_gem("reversible", reversible)
split <- .rc_corda2_split_original(reversible)
stopifnot(
  identical(colnames(split$S), c("REV", "IRR", "REV_CORDA_rev_rxn")),
  split$ub[["REV_CORDA_rev_rxn"]] == 1000
)
engine <- .rc_corda_new_lp_engine(split, "highs", 30)
constrained <- .rc_corda2_constrain_target(
  engine, split, "REV", options
)
engine <- .rc_corda_release_lp_engine(constrained$engine)
stopifnot(
  identical(constrained$opposite, "REV_CORDA_rev_rxn"),
  constrained$upper[["REV_CORDA_rev_rxn"]] == 0,
  constrained$lower[["REV"]] == constrained$upper[["REV"]]
)

network <- make_gem(
  c("A", "B"), c("M", "H", "N"),
  list(
    list("A", "M", 1),
    list("A", "H", -1), list("B", "H", 1),
    list("B", "N", -1)
  )
)
inspect_gem("network", network)
network_split <- .rc_corda2_split_original(network)
result <- .rc_corda_build_three_stage_serial_core(
  split = network_split,
  classes = classes(hc = "H", mc = "M", nc = "N"),
  options = options,
  solver = "highs",
  time_limit = 30
)
stopifnot(
  setequal(result$included, c("M", "H", "N")),
  setequal(result$stage1_associated, c("M", "N")),
  result$HCtoMC["H", "M"] == 1L,
  result$HCtoNC["H", "N"] == 1L,
  identical(result$algorithm,
            "schultzdre_MATLAB_CORDA2_original_semantics"),
  identical(result$stage_update_policy,
            "original_matlab_directional_order"),
  identical(result$parallel_execution_policy,
            "serial_original_persistent_engine")
)

cat("Original MATLAB CORDA2 source-contract checks passed\n")