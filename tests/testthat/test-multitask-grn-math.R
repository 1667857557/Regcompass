test_that("condition-theta design creates one sparse block per condition", {
  x <- cbind(E1 = c(1, 2, 3, 4), E2 = c(5, 6, 7, 8))
  condition <- c("B", "A", "B", "A")
  design <- RegCompassR:::.rc_condition_theta_design_matrix(
    x,
    condition,
    condition_levels = c("A", "B")
  )

  expect_identical(
    colnames(design),
    c("A::E1", "A::E2", "B::E1", "B::E2")
  )
  expect_equal(design[condition == "A", 1:2, drop = FALSE],
               x[condition == "A", , drop = FALSE])
  expect_equal(design[condition == "A", 3:4, drop = FALSE],
               matrix(0, nrow = 2, ncol = 2))
  expect_equal(design[condition == "B", 1:2, drop = FALSE],
               matrix(0, nrow = 2, ncol = 2))
  expect_equal(design[condition == "B", 3:4, drop = FALSE],
               x[condition == "B", , drop = FALSE])
})

test_that("direct theta decoding derives zero-sum deviations", {
  decoded <- RegCompassR:::.rc_decode_condition_theta(
    coefficient = c(
      0.8, 0,
      0.2, -0.6,
      -0.1, 0.3
    ),
    n_edges = 2,
    condition_levels = c("A", "B", "C"),
    edge_ids = c("E1", "E2")
  )

  expect_equal(decoded$theta["A", ], c(E1 = 0.8, E2 = 0))
  expect_equal(decoded$beta, colMeans(decoded$theta))
  expect_equal(colSums(decoded$delta), c(E1 = 0, E2 = 0),
               tolerance = 1e-12)
  expect_equal(decoded$theta,
               sweep(decoded$delta, 2L, decoded$beta, "+"))
  expect_identical(decoded$theta["A", "E2"], 0)
})

test_that("condition balancing gives every condition equal total loss weight", {
  condition <- c(rep("A", 10), rep("B", 25), rep("C", 5))
  weight <- RegCompassR:::.rc_condition_balanced_weights(condition)
  total <- tapply(weight, condition, sum)

  expect_equal(unname(total), rep(total[[1]], 3), tolerance = 1e-12)
  expect_equal(mean(weight), 1, tolerance = 1e-12)
})

test_that("fallback bootstrap resamples full condition sizes by cell", {
  condition <- c(rep("A", 8), rep("B", 11), rep("C", 6))
  set.seed(17)
  index <- RegCompassR:::.rc_condition_stratified_bootstrap_indices(
    condition,
    sample = NULL
  )

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

test_that("sample bootstrap resamples whole donor clusters within condition", {
  condition <- c(rep("A", 6), rep("B", 7))
  sample <- c(
    rep("A1", 2), rep("A2", 1), rep("A3", 3),
    rep("B1", 2), rep("B2", 4), rep("B3", 1)
  )
  set.seed(31)
  index <- RegCompassR:::.rc_sample_cluster_bootstrap_indices(condition, sample)

  expect_setequal(unique(condition[index]), c("A", "B"))
  for (sample_id in unique(sample)) {
    cells <- which(sample == sample_id)
    multiplicity <- vapply(cells, function(cell) {
      sum(index == cell)
    }, integer(1))
    expect_equal(length(unique(multiplicity)), 1L)
  }
  for (level in unique(condition)) {
    sample_ids <- unique(sample[condition == level])
    multiplicity <- vapply(sample_ids, function(sample_id) {
      first_cell <- which(sample == sample_id & condition == level)[[1L]]
      sum(index == first_cell)
    }, integer(1))
    expect_equal(sum(multiplicity), length(sample_ids))
  }
})

test_that("sample bootstrap resolution warns only when it degrades", {
  meta <- data.frame(
    condition = rep(c("A", "B"), each = 4),
    donor = rep(c("d1", "d2", "d3", "d4"), each = 2),
    stringsAsFactors = FALSE
  )
  resolved <- expect_silent(
    RegCompassR:::.rc_resolve_bootstrap_sample(
      meta, sample_col = "donor", condition_col = "condition"
    )
  )
  expect_identical(resolved$resampling_unit, "sample")
  expect_identical(resolved$sample_col, "donor")

  fallback <- expect_warning(
    RegCompassR:::.rc_resolve_bootstrap_sample(
      meta, sample_col = "missing", condition_col = "condition"
    ),
    "metadata column `missing` does not exist"
  )
  expect_identical(fallback$resampling_unit, "cell")
  expect_match(fallback$fallback_reason, "does not exist", fixed = TRUE)
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
    "lasso component"
  )
  expect_error(
    RegCompassR:::.rc_validate_multitask_grn_args(list(alpha = 1)),
    "ridge component"
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

test_that("public canonical workflow exposes sample bootstrap only in Stage 1", {
  expect_true("sample_col" %in% names(formals(rc_run_regcompass)))
  expect_true("sample_col" %in% names(formals(rc_run_regcompass_one_shot)))
  expect_true("sample_col" %in% names(formals(rc_regcompass_step_grn)))
  expect_false("sample_col" %in% names(formals(rc_regcompass_step_metacells)))
})
