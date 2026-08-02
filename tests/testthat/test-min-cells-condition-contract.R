test_that("standard Pando applies min_cells to each cell type", {
  cells_a <- paste0("A_", seq_len(300L))
  cells_b <- paste0("B_", seq_len(299L))
  cells <- c(cells_a, cells_b)
  counts <- Matrix::sparseMatrix(
    i = rep.int(1L, length(cells)),
    j = seq_along(cells),
    x = rep.int(1, length(cells)),
    dims = c(1L, length(cells)),
    dimnames = list("Gene1", cells)
  )
  metadata <- data.frame(
    condition = "single",
    cell_type = c(rep("A", 300L), rep("B", 299L)),
    row.names = cells,
    stringsAsFactors = FALSE
  )
  object <- Seurat::CreateSeuratObject(counts = counts, meta.data = metadata)
  result <- RegCompassR:::.rc_filter_stage1_groups_by_min_cells(
    object = object,
    condition_col = "condition",
    celltype_col = "cell_type",
    cell_type = NULL,
    min_cells = 300L
  )
  expect_identical(result$analysis_mode, "standard_pando")
  expect_identical(result$retained_cell_types, "A")
  expect_equal(ncol(result$object), 300L)
  expect_true(all(result$diagnostics$threshold_scope == "cell_type"))
})

test_that("condition Pando requires min_cells in every condition-celltype stratum", {
  sizes <- c(A_control = 300L, A_treated = 300L,
             B_control = 300L, B_treated = 299L)
  cells <- unlist(Map(
    function(name, n) paste0(name, "_", seq_len(n)),
    names(sizes), sizes,
    USE.NAMES = FALSE
  ))
  cell_type <- c(
    rep("A", sizes[["A_control"]]),
    rep("A", sizes[["A_treated"]]),
    rep("B", sizes[["B_control"]]),
    rep("B", sizes[["B_treated"]])
  )
  condition <- c(
    rep("control", sizes[["A_control"]]),
    rep("treated", sizes[["A_treated"]]),
    rep("control", sizes[["B_control"]]),
    rep("treated", sizes[["B_treated"]])
  )
  counts <- Matrix::sparseMatrix(
    i = rep.int(1L, length(cells)),
    j = seq_along(cells),
    x = rep.int(1, length(cells)),
    dims = c(1L, length(cells)),
    dimnames = list("Gene1", cells)
  )
  metadata <- data.frame(
    condition = condition,
    cell_type = cell_type,
    row.names = cells,
    stringsAsFactors = FALSE
  )
  object <- Seurat::CreateSeuratObject(counts = counts, meta.data = metadata)
  result <- RegCompassR:::.rc_filter_stage1_groups_by_min_cells(
    object = object,
    condition_col = "condition",
    celltype_col = "cell_type",
    cell_type = NULL,
    min_cells = 300L
  )
  expect_identical(result$analysis_mode, "condition_grn")
  expect_identical(result$retained_cell_types, "A")
  expect_equal(ncol(result$object), 600L)
  expect_true(all(result$diagnostics$threshold_scope == "condition_x_cell_type"))
  expect_true(all(result$diagnostics$retained[result$diagnostics$cell_type == "A"]))
  expect_false(any(result$diagnostics$retained[result$diagnostics$cell_type == "B"]))
  expect_equal(
    result$diagnostics$n_cells[
      result$diagnostics$cell_type == "B" &
        result$diagnostics$condition == "treated"
    ],
    299L
  )
})
