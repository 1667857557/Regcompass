test_that("rc_pool_diagnostics reports v0.4 pool diagnostic fields", {
  skip_if_not_installed("Matrix")

  pool_map <- data.frame(
    pool_id = c("pool_1", "pool_1", "pool_2", "pool_2"),
    cell_id = paste0("cell", 1:4),
    low_power_pool = c(FALSE, FALSE, TRUE, TRUE),
    sample_id = c("s1", "s1", "s2", "s2"),
    condition = c("case", "case", "control", "control"),
    cell_type = c("T", "T", "B", "B"),
    seurat_clusters = c("0", "0", "1", "1"),
    stringsAsFactors = FALSE
  )
  rna_counts <- Matrix::Matrix(
    c(1, 0, 2, 0,
      0, 3, 0, 1,
      5, 5, 0, 0),
    nrow = 3,
    ncol = 4,
    sparse = TRUE,
    dimnames = list(c("HK1", "PFKM", "LDHA"), paste0("cell", 1:4))
  )
  out <- rc_pool_diagnostics(
    pool_map,
    rna_counts = rna_counts,
    gpr_genes = c("HK1", "LDHA")
  )

  expect_equal(colnames(out), c(
    "pool_id", "sample_id", "cell_type", "condition", "n_cells",
    "low_power_pool", "single_pool_group_flag", "pool_seed", "state_source",
    "state_resolution", "RNA_depth_mean", "GPR_gene_detection_rate"
  ))
  expect_equal(out$n_cells, c(2L, 2L))
  expect_equal(out$sample_id, c("s1", "s2"))
  expect_true(all(is.finite(out$RNA_depth_mean)))
  expect_true(all(out$GPR_gene_detection_rate >= 0 & out$GPR_gene_detection_rate <= 1))
})

test_that("rc_pool_diagnostics validates matrix cell IDs", {
  skip_if_not_installed("Matrix")

  pool_map <- data.frame(pool_id = "pool_1", cell_id = "cell_missing")
  rna_counts <- Matrix::Matrix(1, nrow = 1, ncol = 1, sparse = TRUE)
  rownames(rna_counts) <- "gene1"
  colnames(rna_counts) <- "cell1"

  expect_error(rc_pool_diagnostics(pool_map, rna_counts = rna_counts), "absent")
})

test_that("feature detection mean handles sparse logical matrices", {
  mat <- Matrix::Matrix(c(1, 0, 0, 2), nrow = 2, sparse = TRUE)
  rownames(mat) <- c("g1", "g2")
  colnames(mat) <- c("c1", "c2")
  expect_equal(rc_pool_diagnostics(
    data.frame(pool_id = "p1", cell_id = colnames(mat), skipped = FALSE),
    rna_counts = mat,
    gpr_genes = rownames(mat)
  )$GPR_gene_detection_rate, 0.5)
})
