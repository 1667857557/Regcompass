test_that("associated() matches Python target and redundancy control flow", {
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_corda2_associated)),
    collapse = "\n"
  )
  expect_match(implementation, ".rc_corda2_penalties", fixed = TRUE)
  expect_match(implementation, "split$ub[[target]] < split$tolerance", fixed = TRUE)
  expect_match(implementation, "iteration < max_iter", fixed = TRUE)
  expect_match(implementation, "options$cost_increase", fixed = TRUE)
  expect_match(implementation, "needed_for_target <- sort(unique", fixed = TRUE)
  expect_match(implementation, "python_serial_target_order", fixed = TRUE)
  expect_match(implementation, "target_parallelism = FALSE", fixed = TRUE)
})

test_that("both direction penalties use forward confidence as in Python", {
  S <- Matrix::Matrix(
    matrix(c(-1, 1), nrow = 2), sparse = TRUE,
    dimnames = list(c("A", "B"), "REV")
  )
  split <- RegCompassR:::.rc_corda_split_model(list(
    S = S, lb = c(REV = -10), ub = c(REV = 10)
  ), tolerance = 1e-7)
  penalty <- RegCompassR:::.rc_corda2_penalties(
    split,
    c("REV::forward" = 3L, "REV::reverse" = -1L),
    penalize_medium = TRUE,
    penalty_factor = 100
  )
  expect_equal(penalty[["REV::forward"]], 0)
  expect_equal(penalty[["REV::reverse"]], 0)

  penalty2 <- RegCompassR:::.rc_corda2_penalties(
    split,
    c("REV::forward" = -1L, "REV::reverse" = 3L),
    penalize_medium = TRUE,
    penalty_factor = 100
  )
  expect_equal(penalty2[["REV::forward"]], 100)
  expect_equal(penalty2[["REV::reverse"]], 100)
})

test_that("reaction confidence reduction equals Python max of both directions", {
  S <- Matrix::Matrix(
    matrix(c(-1, 1, 1, -1), nrow = 2), sparse = TRUE,
    dimnames = list(c("A", "B"), c("R1", "R2"))
  )
  split <- RegCompassR:::.rc_corda_split_model(list(
    S = S,
    lb = c(R1 = -10, R2 = -10),
    ub = c(R1 = 10, R2 = 10)
  ), tolerance = 1e-7)
  directional <- c(
    "R1::forward" = 3L, "R1::reverse" = -1L,
    "R2::forward" = 1L, "R2::reverse" = 2L
  )
  reduced <- RegCompassR:::.rc_corda2_reduce_confidence(
    split, directional
  )
  expect_identical(reduced, c(R1 = 3L, R2 = 2L))
  expect_type(reduced, "integer")
})

test_that("build iteration 2 uses positive-coefficient minimization", {
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_corda2_minimize_medium_targets)),
    collapse = "\n"
  )
  expect_match(implementation, "objective[[target]] <- 1", fixed = TRUE)
  expect_match(
    implementation,
    "answer$objective > options$target_flux",
    fixed = TRUE
  )
  expect_false(grepl("objective[[target]] <- -1", implementation, fixed = TRUE))
})

test_that("direct build function preserves Python serial mutation order", {
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_corda_build_three_stage)),
    collapse = "\n"
  )
  expect_match(implementation, "stage1_targets", fixed = TRUE)
  expect_match(implementation, "stage2_targets", fixed = TRUE)
  expect_match(implementation, "absent_count", fixed = TRUE)
  expect_match(implementation, "split_after_absent$ub", fixed = TRUE)
  expect_match(implementation, ".rc_corda2_minimize_medium_targets", fixed = TRUE)
  expect_match(implementation, "split_stage3$ub[[variable]] <- 0", fixed = TRUE)
  expect_match(implementation, "redundancies = FALSE", fixed = TRUE)
  expect_match(implementation, "python_serial_mutation_order", fixed = TRUE)
})

test_that("persistent HiGHS keeps one model and one-shot fallback", {
  creation <- paste(
    deparse(body(RegCompassR:::.rc_corda_new_lp_engine)), collapse = "\n"
  )
  solve <- paste(
    deparse(body(RegCompassR:::.rc_corda_engine_solve)), collapse = "\n"
  )
  expect_match(creation, "hi_new_solver", fixed = TRUE)
  expect_match(creation, "primal_feasibility_tolerance", fixed = TRUE)
  expect_match(creation, "hi_solver_get_option", fixed = TRUE)
  expect_match(solve, "hi_solver_set_objective", fixed = TRUE)
  expect_match(solve, "hi_solver_set_variable_bounds", fixed = TRUE)
  expect_match(solve, "hi_solver_run", fixed = TRUE)
  expect_match(solve, "one_shot_fallback", fixed = TRUE)
})
