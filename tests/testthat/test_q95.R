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
  expect_equal(out$Q$quantile_used, c(0.90, 0.90))
})
