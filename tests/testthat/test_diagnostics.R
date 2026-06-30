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
  atac_counts <- Matrix::Matrix(
    c(1, 1, 0, 0,
      0, 2, 3, 0),
    nrow = 2,
    ncol = 4,
    sparse = TRUE,
    dimnames = list(c("peak1", "peak2"), paste0("cell", 1:4))
  )

  out <- rc_pool_diagnostics(
    pool_map,
    rna_counts = rna_counts,
    atac_counts = atac_counts,
    state_col = "seurat_clusters",
    metabolic_genes = c("HK1", "PFKM"),
    gpr_genes = c("HK1", "LDHA")
  )

  expect_equal(colnames(out), c(
    "pool_id", "sample_id", "condition", "cell_type", "local_state", "n_cells",
    "low_power_pool", "RNA_depth_mean", "ATAC_depth_mean",
    "metabolic_gene_detection_rate", "GPR_gene_detection_rate"
  ))
  expect_equal(out$n_cells, c(2L, 2L))
  expect_equal(out$sample_id, c("s1", "s2"))
  expect_equal(out$local_state, c("0", "1"))
  expect_true(all(is.finite(out$RNA_depth_mean)))
  expect_true(all(out$metabolic_gene_detection_rate >= 0 & out$metabolic_gene_detection_rate <= 1))
})

test_that("rc_pool_diagnostics validates matrix cell IDs", {
  skip_if_not_installed("Matrix")

  pool_map <- data.frame(pool_id = "pool_1", cell_id = "cell_missing")
  rna_counts <- Matrix::Matrix(1, nrow = 1, ncol = 1, sparse = TRUE)
  rownames(rna_counts) <- "gene1"
  colnames(rna_counts) <- "cell1"

  expect_error(rc_pool_diagnostics(pool_map, rna_counts = rna_counts), "absent")
})
