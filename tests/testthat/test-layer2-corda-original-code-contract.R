test_that("original CORDA2 splits only actively reversible reactions", {
  S <- Matrix::Matrix(
    matrix(c(-1, 1, 1, -1), nrow = 2), sparse = TRUE,
    dimnames = list(c("A", "B"), c("IRR", "REV"))
  )
  split <- RegCompassR:::.rc_corda2_split_original(list(
    S = S,
    lb = c(IRR = 0, REV = -5),
    ub = c(IRR = 10, REV = 7)
  ))
  expect_identical(
    split$direction_table$variable_id,
    c("IRR", "REV", "REV_CORDA_rev_rxn")
  )
  expect_equal(split$ub[["IRR"]], 10)
  expect_equal(split$ub[["REV"]], 7)
  expect_equal(split$ub[["REV_CORDA_rev_rxn"]], 5)
})

test_that("target assessment closes the opposite reversible direction", {
  S <- Matrix::Matrix(
    matrix(c(-1, 1), nrow = 2), sparse = TRUE,
    dimnames = list(c("A", "B"), "REV")
  )
  split <- RegCompassR:::.rc_corda2_split_original(list(
    S = S, lb = c(REV = -10), ub = c(REV = 10)
  ))
  bounds <- RegCompassR:::.rc_corda_target_bounds(
    split, "REV", epsilon = 1
  )
  expect_identical(bounds$opposite_variables, "REV_CORDA_rev_rxn")
  expect_identical(
    bounds$opposite_direction_blocked, "REV_CORDA_rev_rxn"
  )
  expect_equal(bounds$lower[["REV"]], 1)
  expect_equal(bounds$upper[["REV"]], 1)
  expect_equal(bounds$upper[["REV_CORDA_rev_rxn"]], 0)
})

test_that("original CORDA2 parameter defaults match CORDA2.m", {
  options <- RegCompassR:::.rc_layer2_corda_options(list(
    model_completion = "corda2"
  ))
  expect_equal(options$MCxNCthresh, 2)
  expect_equal(options$constraint, 1)
  expect_identical(options$constrainby, "val")
  expect_equal(options$om, 1e4)
  expect_equal(options$ci, 0.01)
  expect_identical(
    options$algorithm,
    "schultzdre_MATLAB_CORDA2_original_semantics"
  )
})

test_that("CORDA2 uses the complete parent without FASTCC", {
  parent_code <- paste(
    deparse(body(RegCompassR:::.rc_corda_parent)), collapse = "\n"
  )
  expect_match(parent_code, "rc_build_full_gem", fixed = TRUE)
  expect_false(grepl("rc_solve_lp", parent_code, fixed = TRUE))
  expect_false(grepl(".rc_fastcc_consistent_reactions", parent_code,
                     fixed = TRUE))
  expect_false(grepl("rc_annotate_reaction_roles", parent_code,
                     fixed = TRUE))
})

test_that("model builder directly invokes the original state machine", {
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_complete_celltype_medium_corda_gem)),
    collapse = "\n"
  )
  expect_match(implementation, ".rc_corda2_split_original", fixed = TRUE)
  expect_match(implementation, ".rc_corda_build_three_stage", fixed = TRUE)
  expect_match(implementation, ".rc_corda2_apply_direction_bounds", fixed = TRUE)
  expect_match(implementation, ".rc_finalize_corda_union_model", fixed = TRUE)
  expect_false(grepl("_base", implementation, fixed = TRUE))
  expect_false(grepl("before_", implementation, fixed = TRUE))
})
