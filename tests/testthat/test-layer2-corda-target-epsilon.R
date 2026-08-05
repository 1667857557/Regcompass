test_that("CORDA2 scoring targets use target flux, not association tolerance", {
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
      corda_options = list(flux_tolerance = 1e-8)
    )
  )
  filtered <- RegCompassR:::.rc_corda2_apply_target_flux(
    model,
    target_flux = 1,
    strict = FALSE,
    cell_type = "A"
  )
  expect_identical(filtered$target_directions$reaction_id, "R1")
  expect_equal(
    filtered$build_params$n_parent_corda2_target_feasible_core_directions,
    2L
  )
  expect_equal(
    filtered$build_params$n_final_corda2_target_feasible_core_directions,
    1L
  )
  expect_identical(
    filtered$closure_diagnostics$completion_status,
    c(
      "corda2_retained_at_target_flux",
      "parent_below_corda2_target_flux",
      "corda2_unresolved_at_target_flux"
    )
  )
  expect_error(
    RegCompassR:::.rc_corda2_apply_target_flux(
      model,
      target_flux = 1,
      strict = TRUE,
      cell_type = "A"
    ),
    "R3:reverse"
  )
})

test_that("CORDA2 keeps association tolerance separate from target flux", {
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_corda2_apply_target_flux)),
    collapse = "\n"
  )
  expect_match(implementation, "vmax >= target_flux", fixed = TRUE)
  expect_match(implementation, "final_vmax >= target_flux", fixed = TRUE)
  expect_match(
    implementation,
    "association_flux_tolerance",
    fixed = TRUE
  )
})

test_that("legacy epsilon helper delegates to CORDA2 target flux", {
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_corda_apply_target_epsilon)),
    collapse = "\n"
  )
  expect_match(implementation, ".rc_corda2_apply_target_flux", fixed = TRUE)
})
