test_that("CORDA-like completion targets every selected HC and MC reaction", {
  implementation <- paste(
    deparse(body(
      RegCompassR:::.rc_complete_celltype_medium_corda_like_gem
    )),
    collapse = "\n"
  )
  expect_match(
    implementation,
    "direction_model, classes$biological",
    fixed = TRUE
  )
  expect_match(
    implementation,
    "corda_completion_target_directions",
    fixed = TRUE
  )
  expect_match(
    implementation,
    "HC/MC directions",
    fixed = TRUE
  )
  expect_match(
    implementation,
    "n_final_feasible_hc_mc_directions",
    fixed = TRUE
  )
})

test_that("CORDA-like scoring targets remain restricted to core reactions", {
  implementation <- paste(
    deparse(body(
      RegCompassR:::.rc_complete_celltype_medium_corda_like_gem
    )),
    collapse = "\n"
  )
  expect_match(
    implementation,
    "scoring_target_directions",
    fixed = TRUE
  )
  expect_match(
    implementation,
    "scoring_parent_feasible_targets",
    fixed = TRUE
  )
  expect_match(
    implementation,
    "final$target_directions <- scoring_parent_feasible_targets",
    fixed = TRUE
  )
})
