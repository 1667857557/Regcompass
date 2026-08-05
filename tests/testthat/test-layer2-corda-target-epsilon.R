test_that("post-build scoring uses the original CORDA2 flux threshold", {
  model <- list(
    closure_diagnostics = data.frame(
      reaction_id = c("R1", "R2", "R3"),
      target_direction = c("forward", "forward", "reverse"),
      feasible = TRUE,
      vmax = c(2, 5e-8, 2),
      final_feasible = TRUE,
      final_vmax = c(2, 5e-8, 5e-8),
      stringsAsFactors = FALSE
    ),
    build_params = list()
  )
  filtered <- RegCompassR:::.rc_corda2_apply_target_flux(
    model,
    flux_threshold = 1e-7,
    strict = TRUE,
    cell_type = "A"
  )
  expect_identical(filtered$target_directions$reaction_id, "R1")
  expect_equal(
    filtered$build_params$n_parent_corda2_feasible_core_directions,
    2L
  )
  expect_equal(
    filtered$build_params$n_final_corda2_feasible_core_directions,
    1L
  )
  expect_identical(
    filtered$closure_diagnostics$completion_status,
    c("corda2_retained", "parent_blocked", "corda2_removed")
  )
  expect_true(filtered$build_params$strict_requested)
  expect_false(filtered$build_params$strict_used_for_reconstruction)
})

test_that("the fixed flux threshold is not a public CORDA2 argument", {
  expect_error(
    RegCompassR:::.rc_layer2_corda_options(list(
      model_completion = "corda2",
      corda2_args = list(fluxThreshold = 1e-6)
    )),
    "Allowed names"
  )
})
