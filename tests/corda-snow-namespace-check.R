suppressPackageStartupMessages({
  library(Matrix)
  library(highs)
  library(BiocParallel)
})

root <- normalizePath(".", mustWork = TRUE)
work <- tempfile("corda2_snow_package_")
pkg <- file.path(work, "Corda2SnowCheck")
lib <- file.path(work, "library")
dir.create(file.path(pkg, "R"), recursive = TRUE)
dir.create(lib, recursive = TRUE)

writeLines(c(
  "Package: Corda2SnowCheck",
  "Title: CORDA2 Snow Namespace Check",
  "Version: 0.0.1",
  "Description: Minimal namespace test for corrected Python CORDA2 workers.",
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
  "    'layer2_corda_output_contract.R'",
  "    'layer2_corda_target_contract.R'",
  "    'layer2_corda_parent_contract.R'",
  "    'layer2_corda2_algorithm.R'",
  "    'run.R'"
), file.path(pkg, "DESCRIPTION"))
writeLines("export(run_corda2_snow_check)", file.path(pkg, "NAMESPACE"))

for (file in c(
  "layer2_corda_evidence.R",
  "layer2_corda_lp.R",
  "layer2_corda_paper_contract.R",
  "layer2_corda_direction_contract.R",
  "layer2_corda_model.R",
  "layer2_corda_output_contract.R",
  "layer2_corda_target_contract.R",
  "layer2_corda_parent_contract.R",
  "layer2_corda2_algorithm.R"
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
  ".rc_layer2_task_bpparam <- function() .TEST_CONTEXT$bpparam",
  ".rc_subset_gem <- function(gem, reactions) gem",
  "rc_prepare_directional_targets <- function(...) data.frame()",
  ".rc_directional_feasibility <- function(...) data.frame()",
  "rc_build_full_gem <- function(gem, ...) gem",
  "rc_export_microcompass <- function(...) invisible(TRUE)"
)
writeLines(support, file.path(pkg, "R", "support.R"))

run <- c(
  "run_corda2_snow_check <- function() {",
  "  metabolites <- c('A', 'B', 'C', 'D', 'E', 'F')",
  "  reactions <- c('SRC_A','M1','SRC_B','M2','H1','SRC_D','N1','M3','M4','SRC_F','M5','N2')",
  "  S <- Matrix::Matrix(0, nrow = length(metabolites), ncol = length(reactions), sparse = TRUE,",
  "                      dimnames = list(metabolites, reactions))",
  "  S['A','SRC_A'] <- 1; S['A','M1'] <- -1; S['C','M1'] <- 1",
  "  S['B','SRC_B'] <- 1; S['B','M2'] <- -1; S['C','M2'] <- 1; S['C','H1'] <- -1",
  "  S['D','SRC_D'] <- 1; S['D','N1'] <- -1; S['E','N1'] <- 1",
  "  S['E','M3'] <- -1; S['E','M4'] <- -1",
  "  S['F','SRC_F'] <- 1; S['F','M5'] <- -1; S['A','N2'] <- -1",
  "  gem <- list(S = S, lb = stats::setNames(rep(0, length(reactions)), reactions),",
  "              ub = stats::setNames(rep(100, length(reactions)), reactions))",
  "  split <- .rc_corda_split_model(gem, tolerance = 1e-8)",
  "  classes <- list(hc = 'H1', mc_module = c('M1','M2','M3','M4','M5'),",
  "                  mc_evidence = character(), mc = c('M1','M2','M3','M4','M5'),",
  "                  nc = c('N1','N2'), ot = c('SRC_A','SRC_B','SRC_D','SRC_F'))",
  "  classes$confidence <- stats::setNames(rep('OT', length(reactions)), reactions)",
  "  classes$confidence[classes$nc] <- 'NC'",
  "  classes$confidence[classes$mc] <- 'MC_module'",
  "  classes$confidence[classes$hc] <- 'HC'",
  "  classes$initial_confidence <- classes$confidence",
  "  options <- .rc_layer2_corda_options(list(model_completion = 'corda2',",
  "    corda2_redundancies = 3L, corda2_penalty_factor = 100,",
  "    corda2_support = 2L, corda2_cost_increase = 1.01,",
  "    corda2_target_flux = 1))",
  "  .TEST_CONTEXT$bpparam <- BiocParallel::SnowParam(",
  "    workers = 2L, type = 'SOCK', progressbar = FALSE,",
  "    exportglobals = TRUE, exportvariables = TRUE)",
  "  on.exit({",
  "    try(BiocParallel::bpstop(.TEST_CONTEXT$bpparam), silent = TRUE)",
  "    .TEST_CONTEXT$bpparam <- FALSE",
  "  }, add = TRUE)",
  "  result <- .rc_corda_build_three_stage(split, classes, options,",
  "                                         solver = 'highs', time_limit = 60)",
  "  stopifnot('N1' %in% result$stage2_promoted_nc,",
  "            all(c('M3','M4','M5') %in% result$stage2_promoted_mc),",
  "            !'N2' %in% result$included,",
  "            identical(result$algorithm,",
  "              'resendislab_python_CORDA2_corrected_redundant_path_assessment'))",
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
suppressPackageStartupMessages(library(Corda2SnowCheck))
stopifnot(isTRUE(run_corda2_snow_check()))
cat("Corrected Python CORDA2 Snow namespace check passed\n")
