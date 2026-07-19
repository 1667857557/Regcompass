test_that("penalty LP rejects incomplete alignment", {
  S <- matrix(c(-1, 1), nrow = 1, dimnames = list("m", c("R1", "R2")))
  expect_error(
    rc_build_abs_penalty_lp(S, c(R1 = 0, R2 = 0), c(R1 = 10, R2 = 10),
                            penalties = 1, target_index = 1, target_min = 1),
    "one value per reaction"
  )
  expect_error(
    rc_build_abs_penalty_lp(S, c(R1 = 0, R2 = 0), c(R1 = 10, R2 = 10),
                            penalties = c(R1 = 1), target_index = 1, target_min = 1),
    "missing reactions"
  )
})

test_that("ranged LP rows are expanded without dropping either bound", {
  A <- Matrix::Matrix(matrix(c(1, 2), nrow = 1), sparse = TRUE)
  out <- .rc_expand_ranged_constraints(A, lhs = 1, rhs = 3)
  expect_equal(nrow(out$A), 2L)
  expect_equal(out$sense, c(">", "<"))
  expect_equal(out$bound, c(1, 3))
})

test_that("suboptimal solver status is not classified as optimal", {
  expect_false(identical(.rc_lp_status("suboptimal"), "optimal"))
})

test_that("multiome penalty uses one expression-derived cost and role overrides", {
  E <- matrix(c(0.5, 0.5), nrow = 2, ncol = 1,
              dimnames = list(c("EX", "R"), "u"))
  roles <- data.frame(
    reaction_id = c("EX", "R"),
    role = c("exchange", "internal"),
    role_source = c("curated", "curated")
  )
  out <- rc_compute_multiome_penalty(E, reaction_roles = roles)
  expect_equal(out$penalty["EX", "u"], 1)
  expect_equal(
    out$penalty["R", "u"],
    1 / (1 + log2(1 + 0.5))
  )
  expect_match(out$evidence_description,
               "integrated into gene support", fixed = TRUE)
})

test_that("blocked full-GEM directions are retained only as diagnostics", {
  S <- matrix(0, nrow = 1, ncol = 2,
              dimnames = list("m", c("blocked", "open")))
  gem <- rc_make_gem(S, lb = c(blocked = 0, open = 0),
                     ub = c(blocked = 0, open = 10))
  directions <- rc_prepare_directional_targets(
    gem, c("blocked", "open"), target_direction = "both"
  )
  allowed <- directions[
    directions$target_direction %in% c("forward", "reverse"),
    , drop = FALSE
  ]
  expect_equal(directions$target_direction[directions$reaction_id == "blocked"],
               "none")
  expect_false("blocked" %in% allowed$reaction_id)
  expect_equal(allowed$target_direction, "forward")
})
