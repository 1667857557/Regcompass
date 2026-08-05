test_that("RegCompass scoring uses fixed Python tflux without changing build", {
  model <- list(
    closure_diagnostics = data.frame(
      reaction_id = c("R1", "R2", "R3"),
      target_direction = c("forward", "forward", "reverse"),
      feasible = TRUE,
      vmax = c(2, 0.5, 2),
      final_feasible = TRUE,
      final_vmax = c(2, 0.5, 0.5),
      stringsAsFactors = FALSE
    ),
    build_params = list(
      corda_options = list(feasibility_tolerance = 1e-7)
    )
  )
  filtered <- RegCompassR:::.rc_corda2_apply_target_flux(
    model,
    target_flux = 1,
    strict = TRUE,
    cell_type = "A"
  )
  expect_identical(filtered$target_directions$reaction_id, "R1")
  expect_equal(
    filtered$build_params$n_parent_corda2_tflux_feasible_core_directions,
    2L
  )
  expect_equal(
    filtered$build_params$n_final_corda2_tflux_feasible_core_directions,
    1L
  )
  expect_identical(
    filtered$closure_diagnostics$completion_status,
    c(
      "corda2_retained_at_tflux",
      "parent_below_corda2_tflux",
      "corda2_unresolved_at_tflux"
    )
  )
  expect_true(filtered$build_params$strict_requested)
  expect_false(filtered$build_params$strict_used_for_reconstruction)
})

test_that("tflux is fixed at one and tolerance remains distinct", {
  expect_error(
    RegCompassR:::.rc_corda2_apply_target_flux(
      list(closure_diagnostics = data.frame()),
      target_flux = 2,
      strict = FALSE,
      cell_type = "A"
    ),
    "fixes `tflux` at 1"
  )
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_corda2_apply_target_flux)),
    collapse = "\n"
  )
  expect_match(implementation, "final_vmax >= 1", fixed = TRUE)
  expect_match(implementation, "strict_used_for_reconstruction", fixed = TRUE)
})
