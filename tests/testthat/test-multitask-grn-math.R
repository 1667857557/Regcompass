test_that("symmetric condition coding has a zero-sum deviation basis", {
  contrast <- RegCompassR:::.rc_multitask_contrast(
    c("treated", "control", "treated", "other")
  )

  expect_equal(rownames(contrast), c("control", "other", "treated"))
  expect_equal(rowSums(contrast), rep(0, 3), tolerance = 1e-12)
  expect_equal(colSums(contrast), rep(0, 3), tolerance = 1e-12)
  expect_equal(unname(contrast), diag(3) - matrix(1 / 3, 3, 3))
})

test_that("decoded condition deviations sum to zero and preserve global mean", {
  condition <- c("A", "B", "C")
  contrast <- RegCompassR:::.rc_multitask_contrast(condition)
  beta <- c(0.5, -0.25)
  gamma <- matrix(
    c(1, 2, -1, 0.5, 0.2, -0.4),
    nrow = 2,
    ncol = 3
  )
  decoded <- RegCompassR:::.rc_decode_multitask_coefficients(
    c(beta, as.numeric(gamma)),
    n_edges = 2,
    contrast = contrast
  )

  expect_equal(colSums(decoded$delta), c(0, 0), tolerance = 1e-12)
  expect_equal(colMeans(decoded$theta), beta, tolerance = 1e-12)
  expect_equal(decoded$theta, sweep(decoded$delta, 2, beta, "+"))
})

test_that("condition balancing gives every condition equal total loss weight", {
  condition <- c(rep("A", 10), rep("B", 25), rep("C", 5))
  weight <- RegCompassR:::.rc_condition_balanced_weights(condition)
  total <- tapply(weight, condition, sum)

  expect_equal(unname(total), rep(total[[1]], 3), tolerance = 1e-12)
  expect_equal(mean(weight), 1, tolerance = 1e-12)
})

test_that("bootstrap resamples full condition sizes with replacement", {
  condition <- c(rep("A", 8), rep("B", 11), rep("C", 6))
  set.seed(17)
  index <- RegCompassR:::.rc_condition_stratified_bootstrap_indices(condition)

  expect_length(index, length(condition))
  expect_equal(
    unname(table(condition[index])),
    unname(table(condition))
  )
  duplicate_within_condition <- vapply(
    split(index, condition[index]),
    function(value) anyDuplicated(value) > 0L,
    logical(1)
  )
  expect_true(any(duplicate_within_condition))
})

test_that("bootstrap data are re-centred within each condition", {
  condition <- c(rep("A", 6), rep("B", 7))
  x <- cbind(x1 = seq_along(condition), x2 = seq_along(condition)^2)
  y <- seq_along(condition) * 0.5
  set.seed(2)
  index <- RegCompassR:::.rc_condition_stratified_bootstrap_indices(condition)
  boot_condition <- condition[index]
  x_centered <- RegCompassR:::.rc_residualize_matrix(x[index, ], boot_condition)
  y_centered <- RegCompassR:::.rc_residualize_vector(y[index], boot_condition)

  for (level in unique(boot_condition)) {
    rows <- boot_condition == level
    expect_equal(
      colMeans(x_centered[rows, , drop = FALSE]), c(0, 0),
      tolerance = 1e-12
    )
    expect_equal(mean(y_centered[rows]), 0, tolerance = 1e-12)
  }
})

test_that("multitask validation requires sparse elastic net and bootstrap quality", {
  expect_error(
    RegCompassR:::.rc_validate_multitask_grn_args(list(alpha = 0)),
    "positive lasso"
  )
  expect_error(
    RegCompassR:::.rc_validate_multitask_grn_args(list(alpha = 1)),
    "non-zero ridge component"
  )
  expect_error(
    RegCompassR:::.rc_validate_multitask_grn_args(list(n_bootstrap = 0)),
    "at least 1"
  )
  expect_error(
    RegCompassR:::.rc_validate_multitask_grn_args(list(
      min_bootstrap_success_fraction = 0
    )),
    "must be in (0, 1]",
    fixed = TRUE
  )
  args <- RegCompassR:::.rc_validate_multitask_grn_args(list(
    alpha = 0.5,
    n_bootstrap = 20,
    min_bootstrap_success_fraction = 0.9
  ))
  expect_equal(args$alpha, 0.5)
  expect_identical(args$n_bootstrap, 20L)
  expect_equal(args$min_bootstrap_success_fraction, 0.9)
  expect_false("n_stability" %in% names(args))
  expect_false("stability_fraction" %in% names(args))
})

test_that("public canonical workflow does not expose a sample column", {
  expect_false("sample_col" %in% names(formals(rc_run_regcompass)))
  expect_false("sample_col" %in% names(formals(rc_regcompass_step_grn)))
  expect_false("sample_col" %in% names(formals(rc_regcompass_step_metacells)))
})
