suppressPackageStartupMessages({
  library(Matrix)
  library(highs)
  library(BiocParallel)
})

root <- normalizePath(".", mustWork = TRUE)
work <- tempfile("corda_snow_package_")
pkg <- file.path(work, "CordaSnowCheck")
lib <- file.path(work, "library")
dir.create(file.path(pkg, "R"), recursive = TRUE)
dir.create(lib, recursive = TRUE)

writeLines(c(
  "Package: CordaSnowCheck",
  "Title: CORDA Snow Namespace Check",
  "Version: 0.0.1",
  "Description: Minimal namespace test for RegCompass CORDA Snow workers.",
  "License: MIT",
  "Encoding: UTF-8",
  "Imports: Matrix, methods, highs, BiocParallel",
  "Collate:",
  "    'support.R'",
  "    'layer2_corda_evidence.R'",
  "    'layer2_corda_lp.R'",
  "    'layer2_corda_paper_contract.R'",
  "    'layer2_corda_direction_contract.R'",
  "    'layer2_corda_model.R'",
  "    'run.R'"
), file.path(pkg, "DESCRIPTION"))
writeLines("export(run_corda_snow_check)", file.path(pkg, "NAMESPACE"))

for (file in c(
  "layer2_corda_evidence.R",
  "layer2_corda_lp.R",
  "layer2_corda_paper_contract.R",
  "layer2_corda_direction_contract.R",
  "layer2_corda_model.R"
)) {
  file.copy(
    file.path(root, "R", file),
    file.path(pkg, "R", file),
    overwrite = TRUE
  )
}

support <- c(
  "`%||%` <- function(x, y) if (is.null(x)) y else x",
  ".rc_as_dgCMatrix <- function(x) methods::as(x, 'dgCMatrix')",
  ".rc_bind_frames_fill <- function(values) {",
  "  values <- values[vapply(values, is.data.frame, logical(1))]",
  "  values <- values[vapply(values, nrow, integer(1)) > 0L]",
  "  if (!length(values)) return(data.frame())",
  "  columns <- unique(unlist(lapply(values, colnames), use.names = FALSE))",
  "  values <- lapply(values, function(value) {",
  "    missing <- setdiff(columns, colnames(value))",
  "    for (name in missing) value[[name]] <- NA",
  "    value[, columns, drop = FALSE]",
  "  })",
  "  answer <- do.call(rbind, values)",
  "  rownames(answer) <- NULL",
  "  answer",
  "}",
  ".rc_lp_status <- function(message = '', code = NA_integer_) {",
  "  text <- tolower(paste(message, collapse = ' '))",
  "  if (grepl('infeasible', text)) return('infeasible')",
  "  if (grepl('unbounded', text)) return('unbounded')",
  "  if (grepl('optimal', text)) return('optimal')",
  "  if (is.finite(code) && as.integer(code) == 0L) return('optimal')",
  "  'error'",
  "}",
  "rc_validate_gem <- function(gem) {",
  "  S <- .rc_as_dgCMatrix(gem$S)",
  "  reactions <- colnames(S)",
  "  lb <- as.numeric(gem$lb[reactions])",
  "  ub <- as.numeric(gem$ub[reactions])",
  "  names(lb) <- names(ub) <- reactions",
  "  list(S = S, lb = lb, ub = ub, reactions = reactions)",
  "}",
  "rc_solve_lp <- function(obj, A, lhs, rhs, lb, ub, solver = 'highs', time_limit = Inf) {",
  "  control <- list(log_to_console = FALSE, threads = 1L, solver = 'simplex')",
  "  if (is.finite(time_limit)) control$time_limit <- time_limit",
  "  answer <- highs::highs_solve(",
  "    L = as.numeric(obj), lower = as.numeric(lb), upper = as.numeric(ub),",
  "    A = A, lhs = as.numeric(lhs), rhs = as.numeric(rhs), maximum = FALSE,",
  "    control = do.call(highs::highs_control, control)",
  "  )",
  "  list(status = .rc_lp_status(answer$status_message, answer$status),",
  "       solution = as.numeric(answer$primal_solution),",
  "       objective = as.numeric(answer$objective_value),",
  "       solver_message = answer$status_message)",
  "}",
  "rc_parallel_config <- function(workers = NULL, backend = 'auto') {",
  "  list(workers = if (is.null(workers)) 2L else as.integer(workers))",
  "}",
  ".rc_with_internal_single_thread <- function(FUN) FUN()",
  "rc_parallel_lapply <- function(X, FUN, BPPARAM = NULL, ...) {",
  "  extra <- list(...) ",
  "  worker_fun <- function(x) do.call(FUN, c(list(x), extra))",
  "  if (identical(BPPARAM, FALSE) || is.null(BPPARAM)) return(lapply(X, worker_fun))",
  "  BiocParallel::bplapply(X, worker_fun, BPPARAM = BPPARAM)",
  "}",
  ".TEST_CONTEXT <- new.env(parent = emptyenv())",
  ".TEST_CONTEXT$bpparam <- FALSE",
  ".rc_layer2_task_bpparam <- function() .TEST_CONTEXT$bpparam"
)
writeLines(support, file.path(pkg, "R", "support.R"))

run <- c(
  "run_corda_snow_check <- function() {",
  "  metabolites <- c('A', 'B', 'X', 'Y')",
  "  reactions <- c('SRC_A', 'NC1', 'HC1', 'NC_SHARED', 'MC1', 'MC2', 'SRC_Y', 'MC3')",
  "  S <- Matrix::Matrix(0, nrow = length(metabolites), ncol = length(reactions), sparse = TRUE,",
  "                      dimnames = list(metabolites, reactions))",
  "  S['A', 'SRC_A'] <- 1; S['A', 'NC1'] <- -1",
  "  S['B', 'NC1'] <- 1; S['B', 'HC1'] <- -1",
  "  S['X', 'NC_SHARED'] <- 1; S['X', 'MC1'] <- -1; S['X', 'MC2'] <- -1",
  "  S['Y', 'SRC_Y'] <- 1; S['Y', 'MC3'] <- -1",
  "  gem <- list(S = S, lb = stats::setNames(rep(0, length(reactions)), reactions),",
  "              ub = stats::setNames(rep(100, length(reactions)), reactions))",
  "  split <- .rc_corda_split_model(gem, tolerance = 1e-8)",
  "  classes <- list(hc = 'HC1', mc_module = c('MC1', 'MC2', 'MC3'),",
  "                  mc_evidence = character(), mc = c('MC1', 'MC2', 'MC3'),",
  "                  nc = c('NC1', 'NC_SHARED'), ot = c('SRC_A', 'SRC_Y'))",
  "  classes$confidence <- stats::setNames(rep('OT', length(reactions)), reactions)",
  "  classes$confidence[classes$nc] <- 'NC'",
  "  classes$confidence[classes$mc] <- 'MC_module'",
  "  classes$confidence[classes$hc] <- 'HC'",
  "  classes$initial_confidence <- classes$confidence",
  "  options <- .rc_layer2_corda_options(list(model_completion = 'corda',",
  "    corda_gamma = 1e5, corda_kappa = 1e-2, corda_epsilon = 1,",
  "    corda_n = 3L, corda_p = 2L, corda_seed = 19L))",
  "  .TEST_CONTEXT$bpparam <- BiocParallel::SnowParam(",
  "    workers = 2L, type = 'SOCK', progressbar = FALSE,",
  "    exportglobals = TRUE, exportvariables = TRUE)",
  "  on.exit({",
  "    try(BiocParallel::bpstop(.TEST_CONTEXT$bpparam), silent = TRUE)",
  "    .TEST_CONTEXT$bpparam <- FALSE",
  "  }, add = TRUE)",
  "  result <- .rc_corda_build_three_stage(split, classes, options,",
  "                                         solver = 'highs', time_limit = 60)",
  "  stopifnot(setequal(result$included, reactions),",
  "            identical(result$stage1_associated, 'NC1'),",
  "            identical(result$stage2_promoted_nc, 'NC_SHARED'),",
  "            setequal(result$stage3_associated_ot, c('SRC_A', 'SRC_Y'))) ",
  "  revS <- Matrix::Matrix(matrix(c(-1, 1), nrow = 2), sparse = TRUE,",
  "    dimnames = list(c('RA', 'RB'), 'REV'))",
  "  rev <- .rc_corda_split_model(list(S = revS, lb = c(REV = -10),",
  "    ub = c(REV = 10)), tolerance = 1e-8)",
  "  bounds <- .rc_corda_target_bounds(rev, 'REV::forward', epsilon = 1)",
  "  stopifnot(identical(bounds$opposite_variables, 'REV::reverse'),",
  "            bounds$upper[['REV::reverse']] == 0)",
  "  TRUE",
  "}"
)
writeLines(run, file.path(pkg, "R", "run.R"))

status <- system2(
  file.path(R.home("bin"), "R"),
  c("CMD", "INSTALL", "--no-byte-compile", "--library", shQuote(lib), shQuote(pkg))
)
stopifnot(identical(status, 0L))
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(CordaSnowCheck))
stopifnot(isTRUE(run_corda_snow_check()))
cat("CORDA Snow namespace check passed\n")
