test_that("target assessment leaves the opposite reversible direction open", {
  S <- Matrix::Matrix(
    matrix(c(-1, 1), nrow = 2), sparse = TRUE,
    dimnames = list(c("A", "B"), "REV")
  )
  split <- RegCompassR:::.rc_corda_split_model(list(
    S = S, lb = c(REV = -10), ub = c(REV = 10)
  ), tolerance = 1e-7)
  bounds <- RegCompassR:::.rc_corda_target_bounds(
    split, "REV::forward", epsilon = 1
  )
  expect_identical(bounds$opposite_variables, "REV::reverse")
  expect_identical(bounds$opposite_direction_blocked, character())
  expect_equal(bounds$lower[["REV::forward"]], 1)
  expect_equal(bounds$upper[["REV::forward"]], 1e6)
  expect_equal(bounds$upper[["REV::reverse"]], 1e6)
})

test_that("constructor and target bound assignments preserve Python order", {
  S <- Matrix::Matrix(
    matrix(c(-1, 1), nrow = 2), sparse = TRUE,
    dimnames = list(c("A", "B"), "SMALL")
  )
  split <- RegCompassR:::.rc_corda_split_model(list(
    S = S, lb = c(SMALL = 0), ub = c(SMALL = 0.5)
  ), tolerance = 1e-7)
  expect_equal(split$ub[["SMALL::forward"]], 1e6)
  split$ub[["SMALL::forward"]] <- 0.5
  expect_error(
    RegCompassR:::.rc_corda_target_bounds(
      split, "SMALL::forward", epsilon = 1
    ),
    "above the current upper bound"
  )

  extreme <- list(
    S = S,
    lb = c(SMALL = -2e6),
    ub = c(SMALL = -1.5e6)
  )
  expect_error(
    RegCompassR:::.rc_corda_split_model(extreme, tolerance = 1e-7),
    "transient lower bound assignment"
  )
})

test_that("CORDA2 uses the complete parent without a feasibility precheck", {
  parent_code <- paste(
    deparse(body(RegCompassR:::.rc_corda_parent)), collapse = "\n"
  )
  expect_match(parent_code, "rc_build_full_gem", fixed = TRUE)
  expect_match(parent_code, "feasibility_precheck", fixed = TRUE)
  expect_false(grepl("rc_solve_lp", parent_code, fixed = TRUE))
  expect_false(grepl(".rc_fastcc_consistent_reactions", parent_code,
                     fixed = TRUE))
  expect_false(grepl("rc_annotate_reaction_roles", parent_code,
                     fixed = TRUE))
})

test_that("FASTCORE parent remains independent of CORDA2", {
  dispatch_code <- paste(
    deparse(body(RegCompassR:::.rc_fastcore_parent)), collapse = "\n"
  )
  expect_false(grepl("corda", dispatch_code, ignore.case = TRUE))
  expect_false(grepl("getOption", dispatch_code, fixed = TRUE))
  expect_false(grepl("_base", dispatch_code, fixed = TRUE))
})

test_that("original model builder directly calls exact CORDA2 state machine", {
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_complete_celltype_medium_corda_gem)),
    collapse = "\n"
  )
  expect_match(implementation, ".rc_corda_parent", fixed = TRUE)
  expect_match(implementation, ".rc_corda_build_three_stage", fixed = TRUE)
  expect_match(implementation, "time_limit = Inf", fixed = TRUE)
  expect_match(implementation, ".rc_corda2_apply_target_flux", fixed = TRUE)
  expect_match(implementation, ".rc_corda_attach_parent_contract", fixed = TRUE)
  expect_match(implementation, ".rc_finalize_corda_union_model", fixed = TRUE)
  expect_match(implementation, "intentional_corrections <- character()", fixed = TRUE)
  expect_false(grepl("_base", implementation, fixed = TRUE))
  expect_false(grepl("before_", implementation, fixed = TRUE))
})

test_that("all aliases select one exact pinned algorithm", {
  for (requested in c("corda2", "corda", "corda_like")) {
    options <- RegCompassR:::.rc_layer2_corda_options(list(
      model_completion = requested
    ))
    expect_identical(options$requested_model_completion, "corda2")
    expect_identical(
      options$algorithm,
      "resendislab_python_CORDA2_c02e06d_exact_semantics"
    )
  }
})
