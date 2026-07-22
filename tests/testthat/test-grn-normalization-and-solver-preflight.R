test_that("normalized assay validation accepts reordered identical cell sets", {
  skip_if_not_installed("Seurat")
  counts <- Matrix::Matrix(
    matrix(c(1, 0, 2, 3, 1, 0), nrow = 2),
    sparse = TRUE,
    dimnames = list(c("G1", "G2"), c("C1", "C2", "C3"))
  )
  object <- Seurat::CreateSeuratObject(counts = counts)
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object[["RNA"]]@data <- object[["RNA"]]@data[, c("C3", "C1", "C2")]

  expect_no_error(
    RegCompassR:::.rc_require_normalized_assay(object, "RNA", "RNA")
  )
  aligned <- RegCompassR:::.rc_align_normalized_assay(object, "RNA", "RNA")
  expect_identical(
    colnames(SeuratObject::GetAssayData(aligned, assay = "RNA", slot = "data")),
    colnames(aligned)
  )
})

test_that("cell-type-shared TF-IDF handles locally absent peaks without warnings", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Signac")
  cells <- c("A1", "A2", "B1", "B2")
  rna <- Matrix::Matrix(
    matrix(c(1, 2, 3, 4, 2, 1, 2, 1), nrow = 2),
    sparse = TRUE,
    dimnames = list(c("G1", "G2"), cells)
  )
  atac <- Matrix::Matrix(
    matrix(c(2, 1, 0, 0, 0, 0, 3, 1), nrow = 2, byrow = TRUE),
    sparse = TRUE,
    dimnames = list(c("chr1-1-10", "chr1-20-30"), cells)
  )
  object <- Seurat::CreateSeuratObject(counts = rna)
  object[["ATAC"]] <- Signac::CreateChromatinAssay(counts = atac, sep = c("-", "-"))
  object$cell_type <- c("A", "A", "B", "B")

  expect_no_warning({
    normalized <- RegCompassR:::.rc_apply_celltype_shared_tfidf(
      object, celltype_col = "cell_type", atac_assay = "ATAC"
    )
  })
  tfidf <- SeuratObject::GetAssayData(normalized, assay = "ATAC", slot = "data")
  expect_true(all(tfidf["chr1-20-30", c("A1", "A2")] == 0))
  expect_true(all(tfidf["chr1-1-10", c("B1", "B2")] == 0))
  expect_equal(
    unname(normalized@misc$regcompass_atac_normalization$n_zero_count_peaks_by_celltype),
    c(1L, 1L)
  )
})

test_that("the default LP solver is a required package with an explicit preflight", {
  description <- read.dcf(file.path(test_path("..", ".."), "DESCRIPTION"))
  imports <- strsplit(description[1, "Imports"], ",")[[1L]]
  imports <- trimws(sub("\\s*\\(.*$", "", imports))
  expect_true("highs" %in% imports)
  expect_no_error(RegCompassR:::.rc_require_lp_solver("highs"))
})

test_that("local Pando installations do not emit the retired remote-metadata warning", {
  source <- paste(
    readLines(test_path("..", "..", "R", "00_meta_module_utils.R"), warn = FALSE),
    collapse = "\n"
  )
  expect_false(grepl("remote metadata are unavailable", source, fixed = TRUE))
  expect_match(source, "local_or_offline_source_api_verified", fixed = TRUE)
})
