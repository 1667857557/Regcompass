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

test_that("Stage 1 TF and peak candidates require five percent detection in any condition", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")

  cells <- paste0("cell", seq_len(80))
  meta <- data.frame(
    condition = rep(c("A", "B"), each = 40),
    cell_type = rep("T", 80),
    row.names = cells,
    stringsAsFactors = FALSE
  )
  rna <- Matrix::Matrix(
    0, nrow = 4, ncol = 80, sparse = TRUE,
    dimnames = list(c("TARGET", "TF5", "TFLOW", "TFB"), cells)
  )
  rna["TF5", c("cell1", "cell2")] <- 1
  rna["TFLOW", "cell1"] <- 1
  rna["TFB", c("cell41", "cell42")] <- 1
  object <- SeuratObject::CreateSeuratObject(counts = rna, meta.data = meta)

  atac <- Matrix::Matrix(
    0, nrow = 3, ncol = 80, sparse = TRUE,
    dimnames = list(c("peak5", "peakLow", "peakB"), cells)
  )
  atac["peak5", c("cell1", "cell2")] <- 1
  atac["peakLow", "cell1"] <- 1
  atac["peakB", c("cell41", "cell42")] <- 1
  object[["ATAC"]] <- SeuratObject::CreateAssayObject(counts = atac)

  motif_tfs <- data.frame(
    motif = c("m1", "m2", "m3"),
    tf = c("TF5", "TFLOW", "TFB"),
    stringsAsFactors = FALSE
  )
  filtered <- RegCompassR:::.rc_stage1_filter_pando_detection_features(
    object = object,
    pando_motif_args = list(motif_tfs = motif_tfs),
    condition_col = "condition",
    rna_assay = "RNA",
    atac_assay = "ATAC",
    cell_type = "T"
  )

  expect_setequal(filtered$pando_motif_args$motif_tfs$tf, c("TF5", "TFB"))
  expect_setequal(rownames(filtered$object[["ATAC"]]), c("peak5", "peakB"))
  expect_true(all(c("TARGET", "TF5", "TFLOW", "TFB") %in%
                  rownames(filtered$object[["RNA"]])))
  expect_equal(filtered$diagnostics$tf_threshold, 0.05)
  expect_equal(filtered$diagnostics$peak_threshold, 0.05)
  expect_equal(filtered$diagnostics$n_candidate_tfs, 3L)
  expect_equal(filtered$diagnostics$n_retained_tfs, 2L)
  expect_equal(filtered$diagnostics$n_candidate_peaks, 3L)
  expect_equal(filtered$diagnostics$n_retained_peaks, 2L)
})

test_that("Stage 1 TF and peak detection uses the same five percent boundary for one condition", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")

  cells <- paste0("cell", seq_len(20))
  meta <- data.frame(
    condition = rep("only", 20), cell_type = rep("T", 20),
    row.names = cells, stringsAsFactors = FALSE
  )
  rna <- Matrix::Matrix(
    0, nrow = 2, ncol = 20, sparse = TRUE,
    dimnames = list(c("TFPASS", "TFZERO"), cells)
  )
  rna["TFPASS", "cell1"] <- 1
  object <- SeuratObject::CreateSeuratObject(counts = rna, meta.data = meta)
  atac <- Matrix::Matrix(
    0, nrow = 2, ncol = 20, sparse = TRUE,
    dimnames = list(c("peakPass", "peakZero"), cells)
  )
  atac["peakPass", "cell1"] <- 1
  object[["ATAC"]] <- SeuratObject::CreateAssayObject(counts = atac)

  filtered <- RegCompassR:::.rc_stage1_filter_pando_detection_features(
    object = object,
    pando_motif_args = list(motif_tfs = data.frame(
      motif = c("m1", "m2"), tf = c("TFPASS", "TFZERO"),
      stringsAsFactors = FALSE
    )),
    condition_col = "condition",
    rna_assay = "RNA",
    atac_assay = "ATAC",
    cell_type = "T"
  )
  expect_identical(filtered$pando_motif_args$motif_tfs$tf, "TFPASS")
  expect_identical(rownames(filtered$object[["ATAC"]]), "peakPass")
})
