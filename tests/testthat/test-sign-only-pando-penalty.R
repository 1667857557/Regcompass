test_that("sign-only penalty direction ignores estimate magnitude", {
  estimate <- c(-20, -0.2, 0, 0.1, 15, NA_real_)
  eligible <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)

  expect_equal(
    .rc_sign_only_edge_direction(estimate, eligible),
    c(-1, -1, 0, 1, 1, 0)
  )

  expect_equal(
    .rc_sign_only_edge_direction(c(0.2, 20), c(TRUE, TRUE)),
    c(1, 1)
  )
  expect_equal(
    .rc_sign_only_edge_direction(c(-0.2, -20), c(TRUE, TRUE)),
    c(-1, -1)
  )
})

test_that("paired-cell predictor activity is invariant to positive rescaling", {
  x <- c(0, 0.25, 0.5, 1, 2, 4)
  base <- .rc_sign_only_predictor_activity(x)
  scaled <- .rc_sign_only_predictor_activity(100 * x)

  expect_equal(base$activity, scaled$activity, tolerance = 1e-12)
  expect_equal(100 * base$scale, scaled$scale, tolerance = 1e-10)
  expect_true(all(base$activity >= 0 & base$activity < 1))
})

test_that("sign-only pair activity preserves exact membership aggregation", {
  cells <- paste0("cell", 1:4)
  rna <- matrix(
    c(1, 2, 3, 4),
    ncol = 1,
    dimnames = list(cells, "TF1")
  )
  atac <- matrix(
    c(1, 0.5, 2, 1),
    ncol = 1,
    dimnames = list(cells, "chr1-1-10")
  )
  membership <- data.frame(
    cell_id = cells,
    metacell_id = c("mc1", "mc1", "mc2", "mc2"),
    stringsAsFactors = FALSE
  )
  edge <- data.frame(
    tf = "TF1",
    region = "chr1-1-10",
    stringsAsFactors = FALSE
  )

  observed <- .rc_sign_only_pair_metacell_activity(
    rna = rna,
    atac = atac,
    edge = edge,
    membership = membership,
    cells = cells
  )

  predictor <- as.numeric(rna[, 1]) * as.numeric(atac[, 1])
  transformed <- .rc_sign_only_predictor_activity(predictor)$activity
  expected <- c(
    mc1 = mean(transformed[1:2]),
    mc2 = mean(transformed[3:4])
  )

  expect_equal(
    as.numeric(observed$activity[1, c("mc1", "mc2")]),
    as.numeric(expected),
    tolerance = 1e-12
  )
  expect_equal(observed$n_requested_pairs, 1L)
  expect_equal(observed$n_mapped_pairs, 1L)
})

test_that("sign-only contribution changes only when direction changes", {
  activity <- c(mc1 = 0.2, mc2 = 0.8)
  positive_small <- .rc_sign_only_edge_direction(0.05, TRUE) * activity
  positive_large <- .rc_sign_only_edge_direction(5, TRUE) * activity
  negative <- .rc_sign_only_edge_direction(-5, TRUE) * activity

  expect_equal(positive_small, positive_large)
  expect_equal(negative, -positive_small)
})
