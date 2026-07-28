test_that("direct theta decoding permits condition-specific exact zeros", {
  decoded <- .rc_decode_condition_theta(
    coefficient = c(2, -1, 0, 3),
    n_edges = 2,
    condition_levels = c("A", "B"),
    edge_ids = c("E1", "E2")
  )

  expect_equal(decoded$theta["A", ], c(E1 = 2, E2 = -1))
  expect_equal(decoded$theta["B", ], c(E1 = 0, E2 = 3))
  expect_equal(decoded$beta, c(E1 = 1, E2 = 1))
  expect_equal(colSums(decoded$delta), c(E1 = 0, E2 = 0))
  expect_identical(decoded$theta["B", "E1"], 0)
})

test_that("condition design assigns predictors only to their own condition", {
  x <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("E1", "E2"))
  )
  design <- .rc_condition_theta_design_matrix(
    x,
    condition = c("A", "B", "A"),
    condition_levels = c("A", "B")
  )

  expect_equal(
    design,
    cbind(
      `A::E1` = c(1, 0, 5),
      `A::E2` = c(2, 0, 6),
      `B::E1` = c(0, 3, 0),
      `B::E2` = c(0, 4, 0)
    )
  )
})

test_that("legacy penalty aliases cannot encode separate deviation shrinkage", {
  expect_error(
    .rc_validate_multitask_grn_args(list(
      global_penalty_factor = 1,
      deviation_penalty_factor = 2
    )),
    "must be equal",
    fixed = TRUE
  )
  args <- .rc_validate_multitask_grn_args(list(
    global_penalty_factor = 1.5,
    deviation_penalty_factor = 1.5
  ))
  expect_equal(args$global_penalty_factor, 1.5)
  expect_equal(args$deviation_penalty_factor, 1.5)
})

test_that("candidate ordering is deterministic and labels remain aligned", {
  edges <- data.frame(
    edge_id = c("E2", "E1"),
    tf_feature_id = c("TF2", "TF1"),
    atac_feature_id = c("P2", "P1"),
    target_feature_id = c("G1", "G1"),
    tf = c("TF2", "TF1"),
    region = c("P2", "P1"),
    target = c("G1", "G1"),
    stringsAsFactors = FALSE
  )
  ordered <- .rc_order_target_edges(edges)
  expect_identical(ordered$edge_id, c("E1", "E2"))
  expect_identical(ordered$tf_feature_id, c("TF1", "TF2"))
  expect_identical(ordered$atac_feature_id, c("P1", "P2"))
})

test_that("direct theta fitter is invariant to candidate row order", {
  skip_if_not_installed("glmnet")
  set.seed(7)
  n <- 120
  cells <- paste0("c", seq_len(n))
  condition <- rep(c("A", "B"), each = n / 2)
  tf1 <- runif(n, 0.2, 2)
  tf2 <- runif(n, 0.2, 2)
  p1 <- runif(n, 0.5, 1.5)
  p2 <- runif(n, 0.5, 1.5)
  x1 <- tf1 * p1
  x2 <- tf2 * p2
  target <- ifelse(condition == "A", 2.5 * x1, -2.0 * x2) +
    rnorm(n, sd = 0.05)

  rna <- Matrix::Matrix(
    rbind(TF1 = tf1, TF2 = tf2, G1 = target),
    sparse = TRUE
  )
  atac <- Matrix::Matrix(rbind(P1 = p1, P2 = p2), sparse = TRUE)
  colnames(rna) <- cells
  colnames(atac) <- cells
  meta <- data.frame(
    condition = condition,
    row.names = cells,
    stringsAsFactors = FALSE
  )
  edges <- data.frame(
    edge_id = c("E1", "E2"),
    tf_feature_id = c("TF1", "TF2"),
    atac_feature_id = c("P1", "P2"),
    target_feature_id = c("G1", "G1"),
    tf = c("TF1", "TF2"),
    region = c("P1", "P2"),
    target = c("G1", "G1"),
    stringsAsFactors = FALSE
  )
  args <- .rc_validate_multitask_grn_args(list(
    lambda_rule = "lambda.min",
    nfolds = 3L,
    n_bootstrap = 5L,
    min_selection_frequency = 0,
    min_sign_stability = 0,
    min_bootstrap_success_fraction = 0.2,
    seed = 11L
  ))

  fit_ordered <- .rc_fit_multitask_target_direct(
    edges, "G1", rna, atac, meta, "condition", args
  )
  fit_shuffled <- .rc_fit_multitask_target_direct(
    edges[c(2, 1), ], "G1", rna, atac, meta, "condition", args
  )

  global_columns <- c("edge_id", "global_estimate", "edge_scale")
  condition_columns <- c(
    "condition", "edge_id", "effective_estimate",
    "selection_frequency", "sign_stability"
  )
  expect_equal(
    fit_ordered$global[, global_columns],
    fit_shuffled$global[, global_columns],
    tolerance = 1e-10
  )
  expect_equal(
    fit_ordered$condition[, condition_columns],
    fit_shuffled$condition[, condition_columns],
    tolerance = 1e-10
  )
  expect_true(all(
    fit_ordered$condition$coefficient_parameterization ==
      "direct_condition_theta"
  ))
})

test_that("obsolete shadow-wrapper bindings are not loaded", {
  namespace <- asNamespace("RegCompassR")
  expect_false(exists(
    ".rc_fit_multitask_target_cv", envir = namespace, inherits = FALSE
  ))
  expect_false(exists(
    ".rc_fit_multitask_target_core", envir = namespace, inherits = FALSE
  ))
  expect_true(exists(
    ".rc_fit_multitask_target_direct", envir = namespace, inherits = FALSE
  ))
})
