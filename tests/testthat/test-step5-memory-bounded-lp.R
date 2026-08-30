test_that("Step 5 builds the canonical Step-2 LP without block intermediates", {
  S <- Matrix::Matrix(
    matrix(c(1, -1, 0, 1, 0, -1), nrow = 2),
    sparse = TRUE,
    dimnames = list(c("m1", "m2"), c("r1", "r2", "r3"))
  )
  lb <- c(r1 = 0, r2 = -2, r3 = 0)
  ub <- c(r1 = 10, r2 = 3, r3 = 7)
  vmax_result <- list(feasible = TRUE, vmax = 4, status = "optimal")

  observed <- RegCompassR:::.rc_compass_step2_prepare(
    S, lb, ub, "r2", vmax_result,
    target_direction = "reverse", omega = 0.95
  )

  n <- ncol(S)
  zero <- Matrix::Matrix(0, nrow = nrow(S), ncol = n, sparse = TRUE)
  mass_balance <- cbind(S, zero)
  positive <- Matrix::Matrix(0, nrow = n, ncol = 2L * n, sparse = TRUE)
  negative <- positive
  positive[cbind(seq_len(n), seq_len(n))] <- 1
  positive[cbind(seq_len(n), n + seq_len(n))] <- -1
  negative[cbind(seq_len(n), seq_len(n))] <- -1
  negative[cbind(seq_len(n), n + seq_len(n))] <- -1
  target <- Matrix::Matrix(0, nrow = 1, ncol = 2L * n, sparse = TRUE)
  target[1, 2] <- -1
  expected_A <- rbind(mass_balance, positive, negative, target)

  expect_true(observed$runnable)
  # Matrix dimnames are not consumed by the LP solver; reaction identity is
  # carried explicitly in `observed$reactions`/`template$reactions`. Compare the
  # numerical constraint system exactly while retaining that name contract.
  expect_equal(
    unname(as.matrix(observed$template$A)),
    unname(as.matrix(expected_A)),
    tolerance = 0
  )
  expect_identical(observed$reactions, colnames(S))
  expect_identical(observed$template$reactions, colnames(S))
  expect_equal(observed$template$lhs,
               c(rep(0, nrow(S)), rep(-Inf, 2L * n), 0.95 * 4))
  expect_equal(observed$template$rhs,
               c(rep(0, nrow(S)), rep(0, 2L * n), Inf))
  expect_equal(observed$template$lb, c(lb, rep(0, n)))
  expect_equal(observed$template$ub, c(ub, pmax(abs(lb), abs(ub))))
})

test_that("Step 5 Vmax workers discard unused primal flux before collection", {
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_build_microcompass_vmax_cache_core)),
    collapse = "\n"
  )
  expect_match(implementation, "flux = numeric()", fixed = TRUE)
  expect_match(implementation, "rc_parallel_lapply", fixed = TRUE)
  expect_match(implementation, ".rc_microcompass_vmax_tasks", fixed = TRUE)
})

test_that("Step 5 memory fix preserves worker-level parallel scheduling", {
  vmax_impl <- paste(
    deparse(body(RegCompassR:::.rc_build_microcompass_vmax_cache_core)),
    collapse = "\n"
  )
  step2_impl <- paste(
    deparse(body(RegCompassR:::.rc_run_celltype_microcompass_engine_reaction_core)),
    collapse = "\n"
  )
  expect_match(vmax_impl, "workers <- .rc_microcompass_worker_count", fixed = TRUE)
  expect_match(vmax_impl, "rc_parallel_lapply", fixed = TRUE)
  expect_match(step2_impl, ".rc_step2_model_batches", fixed = TRUE)
  expect_match(step2_impl, "rc_parallel_lapply", fixed = TRUE)
})

test_that("Step-2 sparse builder no longer materializes redundant constraint blocks", {
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_compass_step2_prepare)),
    collapse = "\n"
  )
  expect_match(implementation, "Matrix::sparseMatrix", fixed = TRUE)
  expect_false(grepl("mass_balance <-", implementation, fixed = TRUE))
  expect_false(grepl("positive <-", implementation, fixed = TRUE))
  expect_false(grepl("negative <-", implementation, fixed = TRUE))
  expect_false(grepl("zero <-", implementation, fixed = TRUE))
})