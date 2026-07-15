solver_for_core_workflow_test <- function() {
  if (requireNamespace("highs", quietly = TRUE)) return("highs")
  if (requireNamespace("Rglpk", quietly = TRUE)) return("glpk")
  if (requireNamespace("gurobi", quietly = TRUE)) return("gurobi")
  NA_character_
}

test_that("LP backend preserves ranged constraints", {
  solver <- solver_for_core_workflow_test()
  skip_if(is.na(solver), "No LP solver available")

  answer <- rc_solve_lp(
    obj = 1,
    A = matrix(1, nrow = 1),
    lhs = 2,
    rhs = 3,
    lb = 0,
    ub = 10,
    solver = solver
  )

  expect_equal(answer$status, "optimal")
  expect_equal(answer$solution, 2, tolerance = 1e-7)
  expect_equal(answer$objective, 2, tolerance = 1e-7)
})

test_that("partial GPR matches remain candidates but not hard core", {
  nodes <- data.frame(
    sample_id = "S1",
    module_id = "M1",
    gene = c("G1", "G3")
  )
  gpr <- data.frame(
    reaction_id = c("R1", "R1", "R2"),
    and_group_id = c(1, 1, 1),
    gene = c("G1", "G2", "G3")
  )

  mapped <- rc_map_meta_module_core_reactions(nodes, gpr)

  expect_false(unique(mapped$is_core[mapped$reaction_id == "R1"]))
  expect_true(unique(mapped$is_partial_candidate[mapped$reaction_id == "R1"]))
  expect_true(unique(mapped$is_core[mapped$reaction_id == "R2"]))
  expect_true(all(c(
    "required_genes", "matched_genes", "missing_genes", "group_complete"
  ) %in% colnames(mapped)))
})
