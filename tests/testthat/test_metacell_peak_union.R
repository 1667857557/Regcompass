test_that("stratum-specific peak matrices are unioned with zero fill", {
  m1 <- Matrix::Matrix(
    matrix(
      c(1, 2, 3, 4),
      nrow = 2,
      dimnames = list(c("peak_a", "peak_b"), c("mc1", "mc2"))
    ),
    sparse = TRUE
  )
  m2 <- Matrix::Matrix(
    matrix(
      c(5, 6, 7, 8),
      nrow = 2,
      dimnames = list(c("peak_b", "peak_c"), c("mc3", "mc4"))
    ),
    sparse = TRUE
  )

  out <- .rc_cbind_sparse_feature_union(list(m1, m2))

  expect_identical(rownames(out), c("peak_a", "peak_b", "peak_c"))
  expect_identical(colnames(out), c("mc1", "mc2", "mc3", "mc4"))
  expect_equal(as.matrix(out), matrix(
    c(
      1, 2, 0,
      3, 4, 0,
      0, 5, 6,
      0, 7, 8
    ),
    nrow = 3,
    dimnames = list(
      c("peak_a", "peak_b", "peak_c"),
      c("mc1", "mc2", "mc3", "mc4")
    )
  ))
})

test_that("peak union rejects duplicated metacell IDs", {
  m1 <- Matrix::Matrix(
    matrix(1, nrow = 1, dimnames = list("peak_a", "mc1")),
    sparse = TRUE
  )
  m2 <- Matrix::Matrix(
    matrix(2, nrow = 1, dimnames = list("peak_b", "mc1")),
    sparse = TRUE
  )

  expect_error(
    .rc_cbind_sparse_feature_union(list(m1, m2)),
    "Duplicated metacell IDs"
  )
})
