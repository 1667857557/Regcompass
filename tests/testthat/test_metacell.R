test_that("metacell metadata and validation enforce unique sample-aware units", {
  membership <- data.frame(
    cell_id = paste0("c", 1:4),
    metacell_id = c("mc1", "mc1", "mc2", "mc2"),
    sample_id = "s1",
    condition = "ctrl",
    cell_type = "T",
    stringsAsFactors = FALSE
  )
  meta <- rc_build_metacell_metadata(membership)
  expect_equal(meta$n_cells, c(2L, 2L))
  counts <- Matrix::Matrix(matrix(1:6, nrow = 3, dimnames = list(paste0("g", 1:3), c("mc1", "mc2"))), sparse = TRUE)
  expect_true(rc_validate_metacell_inputs(counts, meta))
  bad_meta <- rbind(meta, meta[1, ])
  expect_error(rc_validate_metacell_inputs(counts, bad_meta), "unique")
})

test_that("metacell layer1 runs from raw metacell counts", {
  counts <- Matrix::Matrix(
    matrix(c(10, 1, 5, 5, 1, 10), nrow = 3, dimnames = list(c("HK1", "HK2", "PFKM"), c("mc1", "mc2"))),
    sparse = TRUE
  )
  meta <- data.frame(metacell_id = c("mc1", "mc2"), sample_id = c("s1", "s1"), condition = c("ctrl", "ctrl"), cell_type = c("T", "T"), n_cells = c(30L, 25L), stringsAsFactors = FALSE)
  gpr <- data.frame(reaction_id = c("R_HEX", "R_HEX", "R_PFK"), and_group_id = c(1, 2, 1), gene = c("HK1", "HK2", "PFKM"), stringsAsFactors = FALSE)
  out <- rc_run_layer1_from_metacells(gpr, counts, meta, bootstrap = FALSE)
  expect_equal(colnames(out$C_raw), c("mc1", "mc2"))
  expect_equal(out$layer1_unit, "metacell")
  expect_true(all(out$metacell_meta$metacell_id %in% c("mc1", "mc2")))
})

test_that("metacell sample summary reports metacell diagnostics", {
  score <- matrix(c(1, 3, 2, 4), nrow = 2, dimnames = list(c("r1", "r2"), c("mc1", "mc2")))
  meta <- data.frame(metacell_id = c("mc1", "mc2"), sample_id = "s1", condition = "ctrl", cell_type = "T", n_cells = c(20L, 30L), stringsAsFactors = FALSE)
  out <- rc_metacell_sample_summary(score, meta, condition_col = "condition")
  expect_true(all(c("n_metacells_used", "single_metacell_group_flag") %in% colnames(out)))
  expect_equal(unique(out$n_metacells_used), 2L)
  expect_false(unique(out$single_metacell_group_flag))
})
