test_that("rc_sample_aggregate uses sample-condition-celltype medians", {
  score_mat <- matrix(
    c(1, 3, 10, 20,
      2, 4, 30, 50),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("R1", "R2"), paste0("P", 1:4))
  )
  unit_meta <- data.frame(
    pool_id = paste0("P", 1:4),
    sample_id = c("S1", "S1", "S2", "S2"),
    condition = c("ctrl", "ctrl", "stim", "stim"),
    cell_type = c("T", "T", "T", "B"),
    stringsAsFactors = FALSE
  )

  out <- rc_sample_aggregate(score_mat, unit_meta, condition_col = "condition")

  expect_equal(colnames(out), c("S1.ctrl.T", "S2.stim.B", "S2.stim.T"))
  expect_equal(out[, "S1.ctrl.T"], c(R1 = 2, R2 = 3))
  expect_equal(out[, "S2.stim.T"], c(R1 = 10, R2 = 30))
})

test_that("rc_sample_summary reports sample-level flags", {
  score_mat <- matrix(c(1, 3, 2, 4), nrow = 2, dimnames = list(c("R1", "R2"), c("P1", "P2")))
  unit_meta <- data.frame(pool_id = c("P1", "P2"), sample_id = "S1", condition = "ctrl", cell_type = "T", n_cells = c(40, 20), low_power_pool = c(FALSE, TRUE))
  out <- rc_sample_summary(score_mat, unit_meta, condition_col = "condition")
  expect_true(all(c("n_pools_used", "n_cells_used", "low_power_group_flag", "single_pool_group_flag") %in% colnames(out)))
  expect_equal(unique(out$n_pools_used), 2L)
  expect_true(unique(out$low_power_group_flag))
})
