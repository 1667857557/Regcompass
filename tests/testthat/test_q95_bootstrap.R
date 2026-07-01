test_that("rc_q95_bootstrap returns NA for too few finite values", {
  out <- rc_q95_bootstrap(seq_len(10), B = 10)
  expect_true(all(is.na(out)))
  expect_equal(names(out), c("q95", "ci_low", "ci_high", "width"))
})

test_that("rc_q95_bootstrap returns finite CI for enough values", {
  set.seed(1)
  out <- rc_q95_bootstrap(seq_len(30), B = 20)
  expect_true(all(is.finite(out)))
  expect_gte(out[["ci_high"]], out[["ci_low"]])
  expect_gte(out[["width"]], 0)
})

test_that("rc_q95_calibrate includes bootstrap diagnostics", {
  set.seed(1)
  C_raw <- rbind(r1 = seq_len(30), r2 = seq(0, 1, length.out = 30))
  out <- rc_q95_calibrate(C_raw, min_direct = 20, bootstrap = TRUE, B = 20)
  expect_true(all(c("q95_bootstrap", "q95_ci_low", "q95_ci_high", "q95_ci_width") %in% colnames(out$Q)))
  expect_true(all(is.finite(out$Q$q95_ci_width)))
})

test_that("rc_reaction_confidence distinguishes missing detection input", {
  gprs <- list(r1 = list(c("g1", "g2")))
  out <- rc_reaction_confidence(gprs, pool_detection = NULL)
  expect_false(out$detection_available)
  expect_equal(out$missing_gene_fraction, 0)
  expect_true(is.na(out$mean_gpr_detection_rate))
})

test_that("rc_q95_calibrate bootstraps reaction-stratum rows", {
  set.seed(1)
  C_raw <- rbind(r1 = 1:40, r2 = c(rep(NA_real_, 20), 1:20))
  colnames(C_raw) <- paste0("p", 1:40)
  pool_meta <- data.frame(pool_id = colnames(C_raw), cell_type = rep(c("A", "B"), each = 20))
  out <- rc_q95_calibrate(C_raw, pool_meta = pool_meta, stratum_col = "cell_type", bootstrap = TRUE, B = 20)
  expect_equal(nrow(out$Q), 4L)
  expect_true(all(is.na(out$Q$q95_ci_width[out$Q$reaction_id == "r2" & out$Q$stratum == "A"])))
})
