test_that("Stage 1 metabolic targets use an any-condition 20 percent detection union", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")

  cells <- paste0("cell", seq_len(40))
  genes <- c("g20", "G10", "G30", "GOTHER")
  counts <- Matrix::Matrix(
    0, nrow = length(genes), ncol = length(cells), sparse = TRUE,
    dimnames = list(genes, cells)
  )
  counts["g20", c("cell1", "cell2")] <- 1
  counts["G10", c("cell1", "cell11")] <- 1
  counts["G30", c("cell11", "cell12", "cell13")] <- 1
  counts["GOTHER", paste0("cell", 21:25)] <- 1

  meta <- data.frame(
    condition = rep(c("A", "B", "A", "B"), each = 10),
    cell_type = rep(c("T1", "T2"), each = 20),
    row.names = cells,
    stringsAsFactors = FALSE
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = meta)
  gem <- list(metabolic_genes = c("G20", "G10", "G30", "GOTHER"))

  t1 <- subset(object, cells = rownames(meta)[meta$cell_type == "T1"])
  filtered <- RegCompassR:::.rc_stage1_filter_gem_metabolic_genes(
    object = t1, gem = gem, condition_col = "condition",
    rna_assay = "RNA", cell_type = "T1"
  )

  expect_setequal(filtered$metabolic_genes, c("g20", "G30"))
  diagnostics <- attr(filtered, "regcompass_stage1_metabolic_detection")
  expect_equal(diagnostics$threshold, 0.20)
  expect_identical(diagnostics$positive_rule, "RNA counts > 0")
  expect_identical(diagnostics$scope, "broad_cell_type_condition_union")
  expect_equal(unname(diagnostics$condition_detection["g20", "A"]), 0.20)
  expect_equal(unname(diagnostics$condition_detection["G10", ]), c(0.10, 0.10))
  expect_equal(unname(diagnostics$condition_detection["G30", "B"]), 0.30)
  expect_equal(unname(diagnostics$max_detection["GOTHER"]), 0)
})

test_that("Stage 1 metabolic detection reduces to the same 20 percent rule for one condition", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")

  cells <- paste0("cell", seq_len(10))
  counts <- Matrix::Matrix(
    0, nrow = 2, ncol = 10, sparse = TRUE,
    dimnames = list(c("PASS", "FAIL"), cells)
  )
  counts["PASS", c("cell1", "cell2")] <- 1
  counts["FAIL", "cell1"] <- 1
  meta <- data.frame(
    condition = rep("only", 10), cell_type = rep("T", 10),
    row.names = cells, stringsAsFactors = FALSE
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = meta)
  gem <- list(metabolic_genes = c("PASS", "FAIL"))

  filtered <- RegCompassR:::.rc_stage1_filter_gem_metabolic_genes(
    object = object, gem = gem, condition_col = "condition",
    rna_assay = "RNA", cell_type = "T"
  )
  expect_identical(filtered$metabolic_genes, "PASS")
})
