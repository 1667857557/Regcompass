test_that("condition effects are aligned to the shared edge dictionary", {
  condition <- data.frame(
    tf = c("TF1", "TF2"),
    target = c("G1", "G1"),
    region = c("chr1-1-2", "chr1-3-4"),
    term = c("p1:TF1", "p2:TF2"),
    estimate = c(2, -1),
    corr = c(0.4, -0.2),
    stringsAsFactors = FALSE
  )
  universal <- condition[c(2, 1), , drop = FALSE]
  universal$estimate <- c(-0.25, 0.5)
  universal$corr <- c(-0.1, 0.2)
  condition_fit <- data.frame(target = "G1", rsq = 0.8)
  universal_fit <- data.frame(target = "G1", rsq = 0.6)

  observed <- .rc_build_condition_effect_table(
    condition, universal, condition_fit, universal_fit
  )

  expect_equal(observed$condition_estimate, c(2, -1))
  expect_equal(observed$universal_estimate, c(0.5, -0.25))
  expect_equal(observed$condition_effect, c(1.5, -0.75))
  expect_equal(observed$rsq, c(0.8, 0.8))
  expect_equal(observed$universal_rsq, c(0.6, 0.6))
})

test_that("condition and universal networks must share one edge dictionary", {
  condition <- data.frame(
    tf = "TF1", target = "G1", region = "chr1-1-2",
    term = "p1:TF1", estimate = 1, corr = 0.1
  )
  universal <- condition
  universal$region <- "chr1-8-9"
  fit <- data.frame(target = "G1", rsq = 0.5)

  expect_error(
    .rc_build_condition_effect_table(condition, universal, fit, fit),
    "one edge dictionary"
  )
})

test_that("TF-by-ATAC activity uses the Pando interaction predictor", {
  tf <- matrix(c(1, 2, 3, 4), nrow = 2,
               dimnames = list(c("TF1", "TF2"), c("u1", "u2")))
  peak <- matrix(c(5, 6, 7, 8), nrow = 2,
                 dimnames = list(c("P1", "P2"), c("u1", "u2")))

  expect_equal(.rc_tf_peak_interaction(tf, peak), tf * peak)
  expect_error(
    .rc_tf_peak_interaction(tf, peak[, 1, drop = FALSE]),
    "identical dimensions"
  )
})

test_that("regulatory integration remains bounded and zero preserving", {
  rna <- matrix(c(0, 0.5, 1), nrow = 1,
                dimnames = list("g1", c("u1", "u2", "u3")))
  modifier <- matrix(c(1, 1, -1), nrow = 1,
                     dimnames = dimnames(rna))
  observed <- .rc_integrate_regulatory_support(rna, modifier, alpha = 1)

  expect_equal(observed[1, 1], 0)
  expect_equal(observed[1, 3], 1)
  expect_true(observed[1, 2] > 0.5)
  expect_true(all(observed >= 0 & observed <= 1))
})
