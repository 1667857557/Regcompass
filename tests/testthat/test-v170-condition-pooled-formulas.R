test_that("v1.7.0 regulatory integration is bounded and zero preserving", {
  C <- matrix(
    c(0, 0.20, 0.20, 0.80),
    nrow = 4,
    dimnames = list(paste0("g", 1:4), "mc1")
  )
  R <- matrix(
    c(1, 1, -1, 1),
    nrow = 4,
    dimnames = dimnames(C)
  )
  out <- .rc_integrate_regulatory_support_v170(C, R, alpha = 1)

  expect_equal(out["g1", "mc1"], 0)
  expect_gt(out["g2", "mc1"], C["g2", "mc1"])
  expect_lt(out["g3", "mc1"], C["g3", "mc1"])
  expect_gt(out["g4", "mc1"], C["g4", "mc1"])
  expect_true(all(out >= 0 & out <= 1))
})

test_that("v1.7.0 reaction penalty is positive and decreases with expression", {
  E <- matrix(
    c(0, 1, 3, NA_real_),
    nrow = 4,
    dimnames = list(paste0("R", 1:4), "mc1")
  )
  answer <- rc_compute_multiome_penalty(E)
  P <- answer$penalty[, "mc1"]

  expect_equal(P[["R1"]], 1)
  expect_gt(P[["R1"]], P[["R2"]])
  expect_gt(P[["R2"]], P[["R3"]])
  expect_equal(P[["R4"]], 1)
  expect_true(all(is.finite(P) & P > 0))
  expect_identical(answer$penalty_formula, "1 / (1 + log2(1 + E_multiome))")
})

test_that("condition-pooled grouping excludes biological sample", {
  expect_identical(
    .rc_condition_group_cols("condition", "cell_type"),
    c("condition", "cell_type")
  )
  expect_error(
    .rc_condition_group_cols("condition", "condition"),
    "must be distinct"
  )
})
