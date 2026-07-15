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
  out <- rc_q95_calibrate(C_raw, bootstrap = TRUE, B = 20)
  expect_true(all(c("q95_bootstrap", "q95_ci_low", "q95_ci_high", "q95_ci_width") %in% colnames(out$Q)))
  expect_true(all(is.finite(out$Q$q95_ci_width)))
})

test_that("rc_reaction_confidence distinguishes missing detection input", {
  gprs <- list(r1 = list(c("g1", "g2")))
  expect_error(rc_reaction_confidence(gprs, unit_detection = NULL), "Provide `gene_confidence` or `unit_detection`")
})

test_that("rc_q95_calibrate bootstraps reaction-stratum rows", {
  set.seed(1)
  C_raw <- rbind(r1 = 1:40, r2 = c(rep(NA_real_, 20), 1:20))
  colnames(C_raw) <- paste0("p", 1:40)
  unit_meta <- data.frame(pool_id = colnames(C_raw), cell_type = rep(c("A", "B"), each = 20))
  out <- rc_q95_calibrate(C_raw, unit_meta = unit_meta, stratum_col = "cell_type", bootstrap = TRUE, B = 20)
  expect_equal(nrow(out$Q), 4L)
  expect_true(all(is.na(out$Q$q95_ci_width[out$Q$reaction_id == "r2" & out$Q$stratum == "A"])))
})

test_that("Q95 calibration preserves all-NA reactions as NA", {
  C_raw <- rbind(r_all_na = c(NA_real_, NA_real_), r_finite = c(0, 1))
  colnames(C_raw) <- c("p1", "p2")
  out <- rc_q95_calibrate(C_raw, bootstrap = FALSE)
  expect_true(all(is.na(out$C_rel["r_all_na", ])))
  expect_false(any(out$C_rel["r_all_na", ] == 1, na.rm = TRUE))
  expect_true(out$Q$all_missing_reaction_flag[out$Q$reaction_id == "r_all_na"])
})

test_that("Q95 low-n diagnostics are stratum-specific", {
  C_raw <- rbind(r1 = seq_len(40))
  colnames(C_raw) <- paste0("p", seq_len(40))
  unit_meta <- data.frame(pool_id = colnames(C_raw), cell_type = rep(c("A", "B"), each = 20))
  out <- rc_q95_calibrate(C_raw, unit_meta = unit_meta, stratum_col = "cell_type", bootstrap = FALSE)
  expect_equal(out$Q$n_finite, c(20L, 20L))
  expect_equal(out$Q$n_finite_global, c(40L, 40L))
  expect_false(any(out$Q$low_n_flag))
})
