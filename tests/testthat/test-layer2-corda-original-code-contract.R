test_that("CORDA target split is equivalent to leaving j unsplit", {
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
  expect_equal(forward$lower[["REV::reverse"]], 0)
  expect_equal(forward$upper[["REV::reverse"]], 0)

  reverse <- RegCompassR:::.rc_corda_target_bounds(
    split, "REV::reverse", epsilon = 1
  )
  expect_identical(reverse$opposite_variables, "REV::forward")
  expect_equal(reverse$lower[["REV::reverse"]], 1)
  expect_equal(reverse$lower[["REV::forward"]], 0)
  expect_equal(reverse$upper[["REV::forward"]], 0)
})

test_that("CORDA parent does not inherit FASTCC or role pruning", {
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
  expect_false(grepl("parent$lb[forbidden]", parent_code,
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

test_that("CORDA records FASTCORE controls as unused", {
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_complete_celltype_medium_corda_gem)),
    collapse = "\n"
  )
  expect_match(
    implementation,
    "input_fastcore_epsilon_ignored_for_corda",
    fixed = TRUE
  )
  expect_match(implementation, "parent_prepruning", fixed = TRUE)
  expect_match(implementation, "parent_role_blocking", fixed = TRUE)
})
