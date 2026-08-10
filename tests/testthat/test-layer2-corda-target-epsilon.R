test_that("post-build scoring uses one reconstructed-GEM CORDA2 closure", {
  model <- list(
    closure_diagnostics = data.frame(
      reaction_id = c("R1", "R2", "R3"),
      target_direction = c("forward", "forward", "reverse"),
      feasible = c(TRUE, TRUE, FALSE),
      vmax = c(2, 5e-8, NA_real_),
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
  expect_equal(filtered$build_params$n_corda2_tested_core_directions, 3L)
  expect_equal(filtered$build_params$n_corda2_feasible_core_directions, 1L)
  expect_identical(
    filtered$closure_diagnostics$completion_status,
    c("corda2_retained", "corda2_blocked", "corda2_blocked")
  )
  expect_true(filtered$build_params$strict_requested)
  expect_false(filtered$build_params$strict_used_for_reconstruction)
})

test_that("CORDA2 closure solves only the reconstructed GEM", {
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_corda_core_closure_core)),
    collapse = "\n"
  )
  matches <- gregexpr(
    ".rc_corda_closure_directional_feasibility",
    implementation,
    fixed = TRUE
  )[[1L]]
  expect_equal(sum(matches > 0), 1L)
  expect_match(
    implementation,
    "rc_prepare_directional_targets(final, core",
    fixed = TRUE
  )
  expect_false(grepl("parent_diagnostics", implementation, fixed = TRUE))
  expect_false(grepl("final_diagnostics", implementation, fixed = TRUE))
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
