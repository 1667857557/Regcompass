.rc_lp_solver_package <- function(solver) {
  switch(
    as.character(solver),
    highs = "highs",
    glpk = "Rglpk",
    gurobi = "gurobi",
    NA_character_
  )
}

.rc_require_lp_solver <- function(solver) {
  solver <- match.arg(as.character(solver), c("highs", "gurobi", "glpk"))
  package <- .rc_lp_solver_package(solver)
  if (!requireNamespace(package, quietly = TRUE)) {
    installation <- switch(
      solver,
      highs = "install.packages('highs')",
      glpk = "BiocManager::install('Rglpk') or install.packages('Rglpk')",
      gurobi = "install and license the gurobi R package"
    )
    stop(
      "LP solver `", solver, "` is unavailable because package `", package,
      "` is not installed. Run ", installation,
      ". This is a solver installation error, not evidence that the medium-constrained GEM is infeasible.",
      call. = FALSE
    )
  }
  invisible(solver)
}
