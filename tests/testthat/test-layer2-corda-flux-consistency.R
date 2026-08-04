test_that("CORDA implements all three paper stages with barriers", {
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_corda_build_three_stage)),
    collapse = "\n"
  )
  expect_match(implementation, "stage1_hc_dependencies", fixed = TRUE)
  expect_match(implementation, "stage2_mc_nc_support", fixed = TRUE)
  expect_match(
    implementation, "stage2_remaining_mc_feasibility", fixed = TRUE
  )
  expect_match(implementation, "stage3_re_ot_dependencies", fixed = TRUE)
  expect_match(
    implementation, "stage2_NC_supports_at_least_p_MC", fixed = TRUE
  )
  expect_match(
    implementation, "barrier_then_union_order_independent", fixed = TRUE
  )
})

test_that("dependency costs match CORDA confidence semantics", {
  S <- Matrix::Diagonal(4)
  dimnames(S) <- list(paste0("M", 1:4), paste0("R", 1:4))
  gem <- list(
    S = S,
    lb = stats::setNames(rep(0, 4), colnames(S)),
    ub = stats::setNames(rep(10, 4), colnames(S))
  )
  split <- RegCompassR:::.rc_corda_split_model(gem)
  confidence <- c(R1 = "RE", R2 = "MC_module", R3 = "NC", R4 = "OT")
  stage1 <- RegCompassR:::.rc_corda_base_cost(
    split, confidence, "stage1_hc_dependencies", gamma = 1e5
  )
  expect_equal(stage1[["R1::forward"]], 0)
  expect_equal(stage1[["R2::forward"]], 1)
  expect_equal(stage1[["R3::forward"]], 1e5)
  expect_equal(stage1[["R4::forward"]], 0)
  stage2 <- RegCompassR:::.rc_corda_base_cost(
    split, confidence, "stage2_mc_nc_support", gamma = 1e5
  )
  expect_equal(stage2[["R2::forward"]], 0)
  expect_equal(stage2[["R3::forward"]], 1e5)
  stage3 <- RegCompassR:::.rc_corda_base_cost(
    split, confidence, "stage3_re_ot_dependencies", gamma = 1e5
  )
  expect_equal(stage3[["R4::forward"]], 1)
})

test_that("randomized dependency repeats are deterministic by task", {
  first <- RegCompassR:::.rc_corda_noise(
    20, seed = 7L, key = c("stage1", "R1", 2L), kappa = 0.01
  )
  second <- RegCompassR:::.rc_corda_noise(
    20, seed = 7L, key = c("stage1", "R1", 2L), kappa = 0.01
  )
  other <- RegCompassR:::.rc_corda_noise(
    20, seed = 7L, key = c("stage1", "R1", 3L), kappa = 0.01
  )
  expect_identical(first, second)
  expect_false(identical(first, other))
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
    deparse(body(RegCompassR:::.rc_complete_celltype_medium_corda_gem)),
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
