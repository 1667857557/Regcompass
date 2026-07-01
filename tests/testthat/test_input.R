test_that("rc_validate_seurat rejects non-Seurat inputs", {
  expect_error(rc_validate_seurat(list()), "must inherit")
})

test_that("rc_validate_seurat and rc_extract_inputs handle annotated multiome objects", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")

  counts <- Matrix::Matrix(c(1, 0, 2, 3, 0, 1), nrow = 3, sparse = TRUE)
  rownames(counts) <- paste0("gene", seq_len(3))
  colnames(counts) <- paste0("cell", seq_len(2))

  object <- SeuratObject::CreateSeuratObject(counts = counts, assay = "RNA")
  object[["ATAC"]] <- SeuratObject::CreateAssayObject(counts = counts)
  object$sample_id <- c("sample1", "sample2")
  object$cell_type <- c("T cell", "B cell")
  object$condition <- c("case", "control")
  embedding <- matrix(c(1, 2, 3, 4), nrow = 2)
  rownames(embedding) <- colnames(object)
  colnames(embedding) <- c("UMAP_1", "UMAP_2")
  object[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding,
    key = "UMAP_",
    assay = "RNA"
  )

  expect_true(rc_validate_seurat(object, condition_col = "condition", embedding = "umap"))
  extracted <- rc_extract_inputs(object, condition_col = "condition", embedding = "umap")
  expect_identical(dim(extracted$rna), dim(counts))
  expect_identical(dim(extracted$atac), dim(counts))
  expect_true(all(c("sample_id", "cell_type", "condition") %in% colnames(extracted$meta)))
  expect_identical(extracted$embedding, embedding)
})

test_that("rc_validate_seurat reports missing requested metadata", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")

  counts <- Matrix::Matrix(c(1, 0, 2, 3), nrow = 2, sparse = TRUE)
  rownames(counts) <- paste0("gene", seq_len(2))
  colnames(counts) <- paste0("cell", seq_len(2))
  object <- SeuratObject::CreateSeuratObject(counts = counts, assay = "RNA")
  object[["ATAC"]] <- SeuratObject::CreateAssayObject(counts = counts)
  object$sample_id <- c("sample1", "sample2")

  expect_error(rc_validate_seurat(object), "Missing metadata columns: cell_type")
})

test_that("rc_check_metadata reports cell counts and state source", {
  meta <- data.frame(sample_id = c("s1", "s1", "s2"), condition = c("a", "a", "b"), batch = c("x", "y", "x"), cell_type = c("T", "B", "T"), state = c("0", "0", "1"))
  out <- rc_check_metadata(meta, condition_col = "condition", batch_col = "batch", state_col = "state", state_source = "manual")
  expect_true(all(c("cell_counts", "na_counts", "condition_batch", "state_record") %in% names(out)))
  expect_equal(out$state_record$state_source, "manual")
})
