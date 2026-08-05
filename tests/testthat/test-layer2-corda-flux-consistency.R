test_that("CORDA2 implements the Python build stages with barriers", {
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_corda_build_three_stage)),
    collapse = "\n"
  )
  expect_match(
    implementation, "corda2_stage1_high_associations", fixed = TRUE
  )
  expect_match(
    implementation, "corda2_stage2_medium_absent_support", fixed = TRUE
  )
  expect_match(
    implementation, "corda2_stage2_independent_medium_flux", fixed = TRUE
  )
  expect_match(
    implementation, "corda2_stage3_free_completion", fixed = TRUE
  )
  expect_match(
    implementation, "python_build_stage_barriers", fixed = TRUE
  )
  expect_match(
    implementation, "as.integer(absent_count) >= options$support",
    fixed = TRUE
  )
})

test_that("CORDA2 penalties match the Python implementation", {
  confidence <- c(
    H = 3L, M = 2L, L = 1L, U = 0L, N = -1L
  )
  with_medium <- RegCompassR:::.rc_corda2_penalties(
    confidence, penalize_medium = TRUE, penalty_factor = 100
  )
  expect_equal(with_medium[["H"]], 0)
  expect_equal(with_medium[["M"]], 1)
  expect_equal(with_medium[["L"]], 1)
  expect_equal(with_medium[["U"]], 0)
  expect_equal(with_medium[["N"]], 100)
  without_medium <- RegCompassR:::.rc_corda2_penalties(
    confidence, penalize_medium = FALSE, penalty_factor = 100
  )
  expect_equal(without_medium[["M"]], 0)
  expect_equal(without_medium[["L"]], 0)
  expect_equal(without_medium[["N"]], 100)
})

test_that("CORDA2 redundancy search increases only newly used penalties", {
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_corda2_associated_target)),
    collapse = "\n"
  )
  expect_match(implementation, "options$cost_increase", fixed = TRUE)
  expect_match(implementation, "setdiff(candidate, needed)", fixed = TRUE)
  expect_match(implementation, "penalty[weighted_new]", fixed = TRUE)
  expect_match(implementation, "iteration < max_iter", fixed = TRUE)
})

test_that("remaining medium confidence uses maximum flux", {
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_corda2_maximize_targets)),
    collapse = "\n"
  )
  expect_match(
    implementation,
    "objective[[bounds$target_index]] <- -1",
    fixed = TRUE
  )
  expect_match(
    implementation,
    "result$target_flux > options$target_flux",
    fixed = TRUE
  )
  expect_false(grepl(
    "objective[[bounds$target_index]] <- 1",
    implementation,
    fixed = TRUE
  ))
})

test_that("persistent native HiGHS path has one-shot fallback", {
  creation <- paste(
    deparse(body(RegCompassR:::.rc_corda_new_lp_engine)),
    collapse = "\n"
  )
  solve <- paste(
    deparse(body(RegCompassR:::.rc_corda_engine_solve)),
    collapse = "\n"
  )
  expect_match(creation, "hi_new_solver", fixed = TRUE)
  expect_match(creation, "highs_persistent_cpp", fixed = TRUE)
  expect_match(solve, "hi_solver_set_objective", fixed = TRUE)
  expect_match(solve, "hi_solver_set_variable_bounds", fixed = TRUE)
  expect_match(solve, "hi_solver_run", fixed = TRUE)
  expect_match(solve, "basis_reuse", fixed = TRUE)
  expect_match(solve, "one_shot_fallback", fixed = TRUE)
})

test_that("final scoring targets remain restricted to HC core", {
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_complete_celltype_medium_corda_gem_before_corda2)),
    collapse = "\n"
  )
  expect_match(implementation, ".rc_corda_core_closure", fixed = TRUE)
  expect_match(
    implementation, "final$target_directions <- closure$feasible_targets",
    fixed = TRUE
  )
  expect_match(
    implementation, "n_celltype_fastcore_support_reactions = 0L",
    fixed = TRUE
  )
})
