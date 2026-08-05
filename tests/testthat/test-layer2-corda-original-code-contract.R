test_that("CORDA2 closes the opposite reversible target direction", {
  S <- Matrix::Matrix(
    matrix(c(-1, 1), nrow = 2), sparse = TRUE,
    dimnames = list(c("A", "B"), "REV")
  )
  split <- RegCompassR:::.rc_corda_split_model(list(
    S = S,
    lb = c(REV = -10),
    ub = c(REV = 10)
  ))
  forward <- RegCompassR:::.rc_corda_target_bounds(
    split, "REV::forward", epsilon = 1
  )
  expect_identical(forward$opposite_variables, "REV::reverse")
  expect_equal(forward$lower[["REV::forward"]], 1)
  expect_equal(forward$upper[["REV::reverse"]], 0)
  reverse <- RegCompassR:::.rc_corda_target_bounds(
    split, "REV::reverse", epsilon = 1
  )
  expect_identical(reverse$opposite_variables, "REV::forward")
  expect_equal(reverse$lower[["REV::reverse"]], 1)
  expect_equal(reverse$upper[["REV::forward"]], 0)
})

test_that("CORDA2 uses the complete medium-constrained parent", {
  parent_code <- paste(
    deparse(body(RegCompassR:::.rc_corda_parent)), collapse = "\n"
  )
  expect_match(parent_code, "rc_build_full_gem", fixed = TRUE)
  expect_match(parent_code, "corda_parent_prepruning", fixed = TRUE)
  expect_match(parent_code, "corda_parent_role_blocking", fixed = TRUE)
  expect_false(grepl(".rc_fastcc_consistent_reactions", parent_code,
                     fixed = TRUE))
  expect_false(grepl("rc_annotate_reaction_roles", parent_code,
                     fixed = TRUE))
})

test_that("FASTCORE fallback remains the captured original implementation", {
  dispatch_code <- paste(
    deparse(body(RegCompassR:::.rc_fastcore_parent)), collapse = "\n"
  )
  expect_match(
    dispatch_code,
    ".rc_fastcore_parent_before_corda_contract",
    fixed = TRUE
  )
  expect_match(
    dispatch_code,
    "RegCompassR.corda_parent_active",
    fixed = TRUE
  )
  expect_match(dispatch_code, ".rc_corda_parent", fixed = TRUE)
})

test_that("CORDA2 records its pinned Python source", {
  core <- paste(
    deparse(body(
      RegCompassR:::.rc_corda_build_three_stage_before_correction_contract
    )),
    collapse = "\n"
  )
  expect_match(
    core,
    "c02e06d50606bf93f23d8f2e6d6ade0e996ca70e",
    fixed = TRUE
  )
  expect_match(
    core,
    "resendislab_python_CORDA2_corrected_redundant_path_assessment",
    fixed = TRUE
  )
})

test_that("CORDA2 records all intentional corrections", {
  wrapper <- paste(
    deparse(body(RegCompassR:::.rc_corda_build_three_stage)),
    collapse = "\n"
  )
  expect_match(
    wrapper,
    "remaining medium confidence flux is maximized",
    fixed = TRUE
  )
  expect_match(
    wrapper,
    "opposite reversible target direction is blocked",
    fixed = TRUE
  )
  expect_match(
    wrapper,
    "penalties use each directional variable's own current confidence",
    fixed = TRUE
  )
  expect_match(wrapper, "directional_penalty_contract", fixed = TRUE)
})

test_that("CORDA2 has no randomized paper-noise dependency", {
  associated <- paste(
    deparse(body(RegCompassR:::.rc_corda2_associated_target)),
    collapse = "\n"
  )
  expect_false(grepl(".rc_corda_noise", associated, fixed = TRUE))
  expect_match(associated, "options$cost_increase", fixed = TRUE)
  expect_match(associated, "options$penalty_factor", fixed = TRUE)
  options <- RegCompassR:::.rc_layer2_corda_options(list(
    model_completion = "corda2"
  ))
  expect_identical(options$random_noise, FALSE)
  expect_identical(options$deterministic, TRUE)
})

test_that("CORDA2 public aliases all select one algorithm", {
  for (requested in c("corda2", "corda", "corda_like")) {
    options <- RegCompassR:::.rc_layer2_corda_options(list(
      model_completion = requested
    ))
    expect_identical(options$model_completion, "corda2")
    expect_identical(options$requested_model_completion, "corda2")
    expect_identical(
      options$algorithm,
      "resendislab_python_CORDA2_corrected_redundant_path_assessment"
    )
  }
})
