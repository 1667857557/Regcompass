test_that("character SuperCell membership is extracted", {
  skip_if_not_installed("Seurat")
  obj <- Seurat::CreateSeuratObject(matrix(1, nrow = 2, ncol = 3))
  obj@misc$membership <- c("MC1", "MC1", "MC2")
  out <- RegCompassR:::.rc_extract_supercell_membership(
    obj,
    original_cells = paste0("cell", 1:3),
    metacell_ids = c("MC1", "MC2")
  )
  expect_equal(nrow(out), 3)
  expect_setequal(out$metacell_id, c("MC1", "MC2"))
})

test_that("walktrap_clusters SuperCell membership is extracted", {
  skip_if_not_installed("Seurat")
  obj <- Seurat::CreateSeuratObject(matrix(1, nrow = 2, ncol = 3))
  obj@misc$walktrap_clusters <- c(1, 1, 2)
  out <- RegCompassR:::.rc_extract_supercell_membership(
    obj,
    original_cells = paste0("cell", 1:3),
    metacell_ids = c("MC1", "MC2")
  )
  expect_equal(nrow(out), 3)
  expect_equal(out$metacell_id, c("MC1", "MC1", "MC2"))
})
