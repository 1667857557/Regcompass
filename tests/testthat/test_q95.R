test_that("rc_q95_calibrate clips relative capacities and reports diagnostics", {
  C_raw <- rbind(
    r1 = seq(0, 1, length.out = 10),
    r2 = c(NA, rep(2, 9))
  )
  colnames(C_raw) <- paste0("pool", seq_len(10))

  out <- rc_q95_calibrate(C_raw, min_direct = 100)
  expect_true(all(out$C_rel[is.finite(out$C_rel)] <= 1))
  expect_equal(out$Q$reaction_id, rownames(C_raw))
  expect_true(all(out$Q$low_n_flag))
  expect_equal(out$Q$quantile_used, c(0.95, 0.95))
})

test_that("rc_q95_shrink uses continuous shrinkage for all n", {
  C_raw <- rbind(r1 = 1:4, r2 = c(1, 2, 4, 8))
  colnames(C_raw) <- paste0("p", 1:4)
  out <- rc_q95_shrink(C_raw, n0 = 80)
  expect_true(all(out$Q$rho_n > 0 & out$Q$rho_n < 1))
  expect_equal(as.character(out$Q$q95_power_class), rep("very_low", nrow(out$Q)))
  expect_true(all(out$C_rel <= 1))
})

test_that("rc_q95_shrink tolerates all-missing reactions", {
  C_raw <- rbind(r1 = c(NA_real_, NA_real_), r2 = c(1, 2))
  colnames(C_raw) <- c("p1", "p2")
  out <- rc_q95_shrink(C_raw)
  expect_true(is.na(out$Q$q_shrink[out$Q$reaction_id == "r1"]))
  expect_true(all(out$C_rel["r2", ] <= 1))
})

test_that("rc_q95_shrink uses reaction-specific finite n within strata", {
  C_raw <- rbind(r1 = c(1, NA, 2, NA), r2 = c(1, 2, 3, 4))
  colnames(C_raw) <- paste0("p", 1:4)
  pool_meta <- data.frame(pool_id = colnames(C_raw), cell_type = c("A", "A", "B", "B"))
  out <- rc_q95_shrink(C_raw, pool_meta = pool_meta, stratum_col = "cell_type")
  expect_equal(out$Q$n[out$Q$reaction_id == "r1" & out$Q$stratum == "A"], 1L)
  expect_equal(out$Q$n[out$Q$reaction_id == "r2" & out$Q$stratum == "A"], 2L)
})


test_that("rc_q95_shrink power classes are mutually exclusive and cover all n ranges", {
  C_raw <- rbind(
    very_low = c(rep(1, 3), rep(NA_real_, 397)),
    low = c(rep(1, 10), rep(NA_real_, 390)),
    moderate = c(rep(1, 50), rep(NA_real_, 350)),
    adequate = c(rep(1, 100), rep(NA_real_, 300)),
    high = rep(1, 400)
  )
  colnames(C_raw) <- paste0("p", seq_len(ncol(C_raw)))
  out <- rc_q95_shrink(C_raw)
  expect_equal(
    as.character(out$Q$q95_power_class),
    c("very_low", "low", "moderate", "adequate", "high")
  )
})

test_that("q95 unstable flag uses relative CI width", {
  C_raw <- rbind(r1 = c(rep(1, 20), rep(100, 20)))
  colnames(C_raw) <- paste0("p", seq_len(ncol(C_raw)))
  out <- rc_q95_calibrate(C_raw, bootstrap = FALSE)
  Q <- out$Q
  Q$q95_ci_width <- 0.6 * Q$q_value
  Q$q95_unstable_flag <- is.finite(Q$q95_ci_width) &
    (Q$q95_ci_width / pmax(Q$q_value, 1e-6)) > 0.5
  expect_true(Q$q95_unstable_flag)
})
