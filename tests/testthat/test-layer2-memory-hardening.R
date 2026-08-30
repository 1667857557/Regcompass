test_that("Step 5 compacts directional vmax results before collection", {
  value <- list(
    feasible = TRUE,
    vmax = 4.5,
    status = "optimal",
    flux = stats::setNames(seq_len(5000), paste0("R", seq_len(5000)))
  )
  compact <- RegCompassR:::.rc_step5_compact_vmax_result(value)
  expect_true(compact$feasible)
  expect_equal(compact$vmax, 4.5)
  expect_identical(compact$status, "optimal")
  expect_length(compact$flux, 0L)
})

test_that("Step 5 direct sparse LP template is algebraically identical", {
  S <- Matrix::Matrix(
    matrix(
      c(
        1, -1, 0,
        0, 1, -1
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("M1", "M2"), c("R1", "R2", "R3"))
    ),
    sparse = TRUE
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
  target <- Matrix::Matrix(0, nrow = 1L, ncol = 2L * n, sparse = TRUE)
  target[1L, 2L] <- -1
  reference <- rbind(mass_balance, positive, negative, target)

  observed <- RegCompassR:::.rc_step5_abs_flux_constraint_matrix(
    S, target_index = 2L, target_sign = -1
  )
  expect_equal(as.matrix(observed), as.matrix(reference), tolerance = 0)
  expect_equal(dim(observed), dim(reference))
})

test_that("memory-hardened Step 2 keeps the canonical LP result", {
  skip_if_not_installed("highs")
  S <- Matrix::Matrix(
    matrix(
      c(1, -1),
      nrow = 1,
      dimnames = list("M", c("UP", "TARGET"))
    ),
    sparse = TRUE
  )
  lb <- c(UP = 0, TARGET = 0)
  ub <- c(UP = 10, TARGET = 10)
  penalties <- c(UP = 0.25, TARGET = 0.75)
  vmax <- rc_compass_vmax_directional(
    S, lb, ub, "TARGET", direction = "forward", solver = "highs"
  )
  prepared <- RegCompassR:::.rc_compass_step2_prepare(
    S, lb, ub, "TARGET",
    RegCompassR:::.rc_step5_compact_vmax_result(vmax),
    target_direction = "forward", omega = 0.95
  )
  engine <- RegCompassR:::.rc_compass_step2_new_engine(
    prepared$template, "highs"
  )
  on.exit(RegCompassR:::.rc_compass_step2_release_engine(engine), add = TRUE)
  solved <- RegCompassR:::.rc_compass_step2_engine_solve(engine, penalties)
  observed <- RegCompassR:::.rc_compass_step2_result(
    prepared$template, solved$answer
  )
  reference <- rc_compass_two_step_lp_directional(
    S, lb, ub, "TARGET", penalties,
    target_direction = "forward", omega = 0.95, solver = "highs"
  )
  expect_identical(observed$feasible, reference$feasible)
  expect_equal(observed$vmax, reference$vmax, tolerance = 1e-10)
  expect_equal(observed$penalty, reference$penalty, tolerance = 1e-9)
})
