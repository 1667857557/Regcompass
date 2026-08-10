test_that("structural reactions use the COMPASS maximum penalty", {
  expression <- matrix(
    c(NA_real_, 0, 3, 100),
    ncol = 1,
    dimnames = list(
      c("EX", "DM", "SK", "ART"),
      "u1"
    )
  )
  roles <- data.frame(
    reaction_id = rownames(expression),
    role = c("exchange", "demand", "sink", "artificial_support"),
    stringsAsFactors = FALSE
  )

  out <- rc_compute_multiome_penalty(
    expression,
    reaction_roles = roles
  )

  expect_equal(unname(out$penalty[, "u1"]), rep(1, 4))
  expect_identical(
    out$structural_reaction_policy,
    "compass_maximum_expression_penalty_for_all_structural_roles"
  )
  expect_lte(max(out$penalty), 1)
})

test_that("CORDA2 is the default Layer 2 completion", {
  defaults <- RegCompassR:::.rc_layer2_corda_options(list())
  expect_identical(defaults$model_completion, "corda2")
  expect_identical(defaults$requested_model_completion, "corda2")
  expect_identical(
    defaults$algorithm,
    "schultzdre_MATLAB_CORDA2_original_semantics"
  )

  fastcore <- RegCompassR:::.rc_layer2_corda_options(list(
    model_completion = "fastcore"
  ))
  expect_identical(fastcore$model_completion, "fastcore")
  expect_identical(fastcore$algorithm, "add_only_compact_FASTCORE")
})

test_that("CORDA2 finalization restores positive parent lower bounds", {
  parent <- list(
    S = Matrix::Matrix(
      matrix(c(-1, 1), nrow = 2),
      sparse = TRUE,
      dimnames = list(c("A", "B"), "MUST_RUN")
    ),
    lb = c(MUST_RUN = 2),
    ub = c(MUST_RUN = 10)
  )
  split <- RegCompassR:::.rc_corda2_split_original(parent)
  final <- RegCompassR:::.rc_corda2_apply_direction_bounds(
    parent = parent,
    included_variables = "MUST_RUN",
    split = split
  )

  expect_equal(final$lb[["MUST_RUN"]], 2)
  expect_equal(final$ub[["MUST_RUN"]], 10)
})

test_that("either retained CORDA2 direction restores reversible parent bounds", {
  parent <- list(
    S = Matrix::Matrix(
      matrix(c(-1, 1), nrow = 2),
      sparse = TRUE,
      dimnames = list(c("A", "B"), "REV")
    ),
    lb = c(REV = -10),
    ub = c(REV = 8)
  )
  split <- RegCompassR:::.rc_corda2_split_original(parent)

  forward_only <- RegCompassR:::.rc_corda2_apply_direction_bounds(
    parent = parent,
    included_variables = "REV",
    split = split
  )
  reverse_only <- RegCompassR:::.rc_corda2_apply_direction_bounds(
    parent = parent,
    included_variables = "REV_CORDA_rev_rxn",
    split = split
  )

  expect_equal(forward_only$lb[["REV"]], -10)
  expect_equal(forward_only$ub[["REV"]], 8)
  expect_equal(reverse_only$lb[["REV"]], -10)
  expect_equal(reverse_only$ub[["REV"]], 8)
})
