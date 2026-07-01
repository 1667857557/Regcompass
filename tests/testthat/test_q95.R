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
  expect_true(all(out$Q$q95_very_low_power))
  expect_true(all(out$C_rel <= 1))
})

test_that("rc_q95_shrink tolerates all-missing reactions", {
  C_raw <- rbind(r1 = c(NA_real_, NA_real_), r2 = c(1, 2))
  colnames(C_raw) <- c("p1", "p2")
  out <- rc_q95_shrink(C_raw)
  expect_true(is.na(out$Q$q_shrink[out$Q$reaction_id == "r1"]))
  expect_true(all(out$C_rel["r2", ] <= 1))
})
