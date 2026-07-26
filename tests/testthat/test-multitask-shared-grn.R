test_that("condition weights give every task equal total loss", {
  condition <- c(rep("A", 2), rep("B", 5), rep("C", 3))
  weight <- RegCompassR:::.rc_mt_condition_weights(condition)
  totals <- tapply(weight, condition, sum)
  expect_equal(unname(totals), rep(totals[[1L]], length(totals)), tolerance = 1e-12)
  expect_equal(mean(weight), 1, tolerance = 1e-12)
})

test_that("canonical global and deviations are reference free and zero sum", {
  coefficients <- c(
    "G::E1" = 0.4,
    "G::E2" = -0.2,
    "D::A::E1" = 0.6,
    "D::A::E2" = 0.1,
    "D::B::E1" = -0.2,
    "D::B::E2" = -0.3
  )
  fit_ab <- RegCompassR:::.rc_mt_extract_effective(
    coefficients, c("E1", "E2"), c("A", "B")
  )
  fit_ba <- RegCompassR:::.rc_mt_extract_effective(
    coefficients, c("E1", "E2"), c("B", "A")
  )

  expect_equal(fit_ab$effective[, c("A", "B")],
               fit_ba$effective[, c("A", "B")])
  expect_equal(fit_ab$global, rowMeans(fit_ab$effective))
  expect_equal(rowSums(fit_ab$deviation), c(E1 = 0, E2 = 0), tolerance = 1e-12)
  expect_equal(
    sweep(fit_ab$deviation, 1L, fit_ab$global, "+"),
    fit_ab$effective
  )
})

test_that("balanced interaction scale averages within-condition variance", {
  X <- Matrix::Matrix(
    cbind(
      c(0, 2, 10, 14),
      c(1, 1, 2, 2)
    ),
    sparse = TRUE
  )
  value <- RegCompassR:::.rc_mt_balanced_interaction_scale(
    X, c("A", "A", "B", "B"), min_scale = 0.25
  )
  expected_first <- sqrt((1 + 4) / 2)
  expect_equal(value$raw[[1L]], expected_first, tolerance = 1e-12)
  expect_equal(value$scale[[2L]], 0.25, tolerance = 1e-12)
})

test_that("shared projection uses coefficient stability TF reference and scale", {
  edges <- data.frame(
    effective_estimate = c(2, -1),
    stability_weight = c(0.5, 0.8),
    tf_reference = c(3, 2),
    interaction_scale = c(2, 4)
  )
  value <- RegCompassR:::.rc_regulatory_edge_projection_weight(edges)
  expect_equal(value, c(1.5, -0.4), tolerance = 1e-12)
})

test_that("sign flips require active stable opposite directions", {
  edges <- data.frame(
    cell_type = rep("T", 4),
    edge_id = c("E1", "E1", "E2", "E2"),
    effective_estimate = c(1, -2, 1, -2),
    active_edge = c(TRUE, TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  value <- RegCompassR:::.rc_mt_mark_sign_flips(
    edges, "cell_type", tolerance = 1e-8
  )
  expect_true(all(value$sign_flip_flag[value$edge_id == "E1"]))
  expect_false(any(value$sign_flip_flag[value$edge_id == "E2"]))
  expect_equal(
    value$effective_direction,
    c("positive", "negative", "positive", "negative")
  )
})

test_that("zero-variance targets return a schema-compatible skipped fit", {
  edge <- data.frame(
    edge_id = "E1",
    tf = "TF1",
    region = "P1",
    target = "G1",
    tf_feature_id = "TF1",
    atac_feature_id = "P1",
    target_feature_id = "G1",
    interaction_scale = 1,
    raw_interaction_scale = 1,
    tf_reference = 1,
    estimable = TRUE,
    stringsAsFactors = FALSE
  )
  target_design <- list(
    edge_metadata = edge,
    condition_levels = c("A", "B"),
    estimable = TRUE
  )
  args <- RegCompassR:::.rc_mt_validate_args(list(n_stability = 0L))
  fit <- RegCompassR:::.rc_mt_fit_target(
    y = rep(1, 8),
    target_design = target_design,
    condition = rep(c("A", "B"), each = 4),
    args = args,
    target = "G1"
  )
  expect_identical(fit$fit_status, "skipped_zero_target_variance")
  expect_equal(dim(fit$effects$effective), c(1L, 2L))
  expect_equal(fit$effects$effective, matrix(
    0, nrow = 1, ncol = 2,
    dimnames = list("E1", c("A", "B"))
  ))
})

test_that("condition target table carries positive and negative active edges", {
  significant <- data.frame(
    group_id = c("A::T", "A::T", "B::T"),
    condition = c("A", "A", "B"),
    cell_type = c("T", "T", "T"),
    target = c("G1", "G1", "G2"),
    effective_estimate = c(1, -0.5, -1),
    stringsAsFactors = FALSE
  )
  value <- RegCompassR:::.rc_mt_condition_target_table(
    significant, "condition", "cell_type"
  )
  g1 <- value[value$target == "G1", , drop = FALSE]
  expect_equal(g1$n_active_edges, 2L)
  expect_equal(g1$n_positive_edges, 1L)
  expect_equal(g1$n_negative_edges, 1L)
})

test_that("RegCompassR release version is 1.8.8", {
  expect_identical(as.character(utils::packageVersion("RegCompassR")), "1.8.8")
})
