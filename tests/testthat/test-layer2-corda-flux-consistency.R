test_that("Step 1 costs match original CORDA2.m", {
  S <- Matrix::Matrix(
    matrix(c(-1, 1, 1, -1, -1, 1), nrow = 2), sparse = TRUE,
    dimnames = list(c("A", "B"), c("H", "M", "N"))
  )
  split <- RegCompassR:::.rc_corda2_split_original(list(
    S = S,
    lb = c(H = 0, M = 0, N = 0),
    ub = c(H = 10, M = 10, N = 10)
  ))
  directional_class <- c(H = "HC", M = "MC", N = "NC")
  options <- RegCompassR:::.rc_layer2_corda_options(list(
    model_completion = "corda2"
  ))
  cost <- RegCompassR:::.rc_corda2_stage_cost(
    split, directional_class, options,
    penalized_class = "stage1",
    baseline = options$baseline_cost
  )
  expect_equal(cost[["H"]], 1e-3)
  expect_equal(cost[["M"]], sqrt(1e4))
  expect_equal(cost[["N"]], 1e4)
})

test_that("dependency assessment uses original target and cost updates", {
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_corda2_dependency_assessment)),
    collapse = "\n"
  )
  expect_match(implementation, ".rc_corda2_constrain_target", fixed = TRUE)
  expect_match(implementation, "options$flux_threshold", fixed = TRUE)
  expect_match(implementation, "penalty[newly_used]", fixed = TRUE)
  expect_match(implementation, "1 + options$ci", fixed = TRUE)
})

test_that("build follows original Step 1, Step 2.1, Step 2.2 and Step 3", {
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_corda_build_three_stage)),
    collapse = "\n"
  )
  expect_match(implementation, "HCtoMC", fixed = TRUE)
  expect_match(implementation, "HCtoNC", fixed = TRUE)
  expect_match(implementation, "MCxNC", fixed = TRUE)
  expect_match(implementation, "options$MCxNCthresh", fixed = TRUE)
  expect_match(implementation, "split_step22$ub[nc] <- 0", fixed = TRUE)
  expect_match(implementation, "corda2_step3_HC_OT_dependencies", fixed = TRUE)
  expect_match(
    implementation,
    "original_matlab_directional_order",
    fixed = TRUE
  )
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
  expect_match(solve, "hi_solver_set_objective", fixed = TRUE)
  expect_match(solve, "hi_solver_set_variable_bounds", fixed = TRUE)
  expect_match(solve, "one_shot_fallback", fixed = TRUE)
})
