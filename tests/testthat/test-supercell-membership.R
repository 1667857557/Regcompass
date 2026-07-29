test_that("canonical SuperCell membership_table is extracted", {
  skip_if_not_installed("Seurat")
  counts <- matrix(
    1, nrow = 2, ncol = 3,
    dimnames = list(paste0("gene", 1:2), paste0("cell", 1:3))
  )
  obj <- Seurat::CreateSeuratObject(counts)
  obj@misc$membership_table <- data.frame(
    cell_id = paste0("cell", 1:3),
    metacell_id = c("MC1", "MC1", "MC2"),
    stringsAsFactors = FALSE
  )
  out <- RegCompassR:::.rc_extract_supercell_membership(
    obj,
    original_cells = paste0("cell", 1:3),
    metacell_ids = c("MC1", "MC2")
  )
  expect_equal(nrow(out), 3)
  expect_setequal(out$metacell_id, c("MC1", "MC2"))
})

test_that("legacy inferred SuperCell membership is rejected", {
  skip_if_not_installed("Seurat")
  counts <- matrix(
    1, nrow = 2, ncol = 3,
    dimnames = list(paste0("gene", 1:2), paste0("cell", 1:3))
  )
  obj <- Seurat::CreateSeuratObject(counts)
  obj@misc$walktrap_clusters <- c(1, 1, 2)
  expect_error(
    RegCompassR:::.rc_extract_supercell_membership(
      obj,
      original_cells = paste0("cell", 1:3),
      metacell_ids = c("MC1", "MC2")
    ),
    "misc\\$membership_table"
  )
})

test_that("metacell metadata may mix biological samples but not biology strata", {
  membership <- data.frame(
    cell_id = c("c1", "c2"),
    metacell_id = c("MC1", "MC1"),
    sample_id = c("S1", "S2"),
    condition = c("Control", "Control"),
    cell_type = c("T", "T"),
    stringsAsFactors = FALSE
  )
  out <- rc_build_metacell_metadata(membership)
  expect_equal(out$n_cells, 2L)
  expect_error(
    rc_build_metacell_metadata(transform(
      membership, cell_type = c("T", "B")
    )),
    "cell_type"
  )
})
