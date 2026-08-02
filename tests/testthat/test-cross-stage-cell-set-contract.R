test_that("analysis cell set excludes undersized strata and non-estimable types", {
  sizes <- c(
    A_control = 300L,
    A_treated = 299L,
    A_recovery = 300L,
    B_control = 300L,
    B_treated = 299L
  )
  cells <- unlist(Map(
    function(name, n) paste0(name, "_", seq_len(n)),
    names(sizes), sizes,
    USE.NAMES = FALSE
  ))
  metadata <- data.frame(
    condition = c(
      rep("control", sizes[["A_control"]]),
      rep("treated", sizes[["A_treated"]]),
      rep("recovery", sizes[["A_recovery"]]),
      rep("control", sizes[["B_control"]]),
      rep("treated", sizes[["B_treated"]])
    ),
    cell_type = c(
      rep("A", sizes[["A_control"]]),
      rep("A", sizes[["A_treated"]]),
      rep("A", sizes[["A_recovery"]]),
      rep("B", sizes[["B_control"]]),
      rep("B", sizes[["B_treated"]])
    ),
    row.names = cells,
    stringsAsFactors = FALSE
  )
  counts <- Matrix::sparseMatrix(
    i = rep.int(1L, length(cells)),
    j = seq_along(cells),
    x = rep.int(1, length(cells)),
    dims = c(1L, length(cells)),
    dimnames = list("Gene1", cells)
  )
  object <- Seurat::CreateSeuratObject(counts = counts, meta.data = metadata)
  result <- RegCompassR:::.rc_build_stage_analysis_cell_set(
    object = object,
    condition_col = "condition",
    celltype_col = "cell_type",
    cell_type = NULL,
    pando_args = list(min_cells = 300L)
  )

  expect_identical(result$retained_cell_types, "A")
  expect_equal(ncol(result$object), 600L)
  expect_true(all(result$object$cell_type == "A"))
  expect_setequal(unique(result$object$condition), c("control", "recovery"))
  expect_false(any(grepl("A_treated|B_", result$retained_cells)))
  expect_identical(result$skipped_condition_cell_types, "B")
})

test_that("Stage 2 reproduces the exact ordered Stage 1 cell IDs", {
  cells <- c("c1", "c2", "c3")
  counts <- Matrix::sparseMatrix(
    i = rep.int(1L, 3L), j = seq_len(3L), x = rep.int(1, 3L),
    dims = c(1L, 3L), dimnames = list("Gene1", cells)
  )
  object <- Seurat::CreateSeuratObject(counts = counts)
  contract <- list(
    source = "stage1_min_cells_before_normalization",
    retained_cells = c("c3", "c1"),
    retained_cell_types = "A",
    min_cells = 300L
  )
  filtered <- RegCompassR:::.rc_subset_to_stage1_cell_set(object, contract)
  expect_identical(colnames(filtered), c("c3", "c1"))
  expect_equal(
    filtered@misc$regcompass_cross_stage_cell_set$n_cells,
    2L
  )

  missing <- contract
  missing$retained_cells <- c("c3", "absent")
  expect_error(
    RegCompassR:::.rc_subset_to_stage1_cell_set(object, missing),
    "missing 1 cell"
  )
})

test_that("public workflow passes the Stage 1 contract to Stage 2", {
  expect_true("grn" %in% names(formals(rc_regcompass_step_metacells)))
  workflow <- paste(deparse(body(rc_run_regcompass)), collapse = "\n")
  expect_match(workflow, "grn = step1", fixed = TRUE)
  stage2 <- paste(deparse(body(rc_regcompass_step_metacells)), collapse = "\n")
  expect_match(stage2, ".rc_subset_to_stage1_cell_set", fixed = TRUE)
  expect_match(stage2, "stage1_exact_cell_ids_v1", fixed = TRUE)
  expect_match(
    stage2,
    '"requested_condition_col" %in% names(grn$params)',
    fixed = TRUE
  )
})
