test_that("post-CORDA2 targets come directly from final GEM bounds", {
  S <- Matrix::Matrix(
    matrix(c(-1, 1, -1, 1), nrow = 2),
    sparse = TRUE,
    dimnames = list(c("A", "B"), c("REV", "FWD"))
  )
  model <- list(
    S = S,
    lb = c(REV = -10, FWD = 0),
    ub = c(REV = 10, FWD = 5),
    build_params = list()
  )

  prepared <- RegCompassR:::.rc_corda2_prepare_scoring_targets(
    model = model,
    core_reactions = c("REV", "FWD"),
    target_direction = "both",
    strict = TRUE,
    cell_type = "A"
  )

  expect_identical(
    prepared$target_directions$reaction_id,
    c("REV", "REV", "FWD")
  )
  expect_identical(
    prepared$target_directions$target_direction,
    c("forward", "reverse", "forward")
  )
  expect_identical(prepared$required_core_reactions, c("REV", "FWD"))
  expect_identical(prepared$target_status, "ready_for_compass_vmax")
  expect_equal(nrow(prepared$closure_diagnostics), 0L)
  expect_false(prepared$build_params$post_reconstruction_closure_lp)
  expect_equal(prepared$build_params$n_required_core_reactions, 2L)
  expect_equal(prepared$build_params$n_retained_core_reactions, 2L)
  expect_equal(prepared$build_params$n_core_scoring_directions, 3L)
})

test_that("CORDA2 builder has no post-reconstruction closure LP", {
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_complete_celltype_medium_corda_gem_core)),
    collapse = "\n"
  )
  expect_match(
    implementation,
    ".rc_corda2_prepare_scoring_targets",
    fixed = TRUE
  )
  expect_match(implementation, "post_build_lp_solves=0", fixed = TRUE)
  expect_false(grepl(".rc_corda_core_closure", implementation, fixed = TRUE))
  expect_false(grepl(
    ".rc_corda_closure_directional_feasibility",
    implementation,
    fixed = TRUE
  ))
  expect_false(grepl(
    ".rc_directional_feasibility_core",
    implementation,
    fixed = TRUE
  ))
})

test_that("missing core reactions are a hard final-model error", {
  model <- list(
    S = Matrix::Matrix(
      matrix(c(-1, 1), nrow = 2),
      sparse = TRUE,
      dimnames = list(c("A", "B"), "R1")
    ),
    lb = c(R1 = 0),
    ub = c(R1 = 10),
    build_params = list()
  )
  expect_error(
    RegCompassR:::.rc_corda2_prepare_scoring_targets(
      model,
      core_reactions = c("R1", "R2"),
      target_direction = "both",
      strict = FALSE,
      cell_type = "A"
    ),
    "immutable structural backbone"
  )
})

test_that("directional vmax feasibility belongs to microCOMPASS scoring", {
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_build_microcompass_vmax_cache_core)),
    collapse = "\n"
  )
  expect_match(implementation, "rc_compass_vmax_directional", fixed = TRUE)
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
