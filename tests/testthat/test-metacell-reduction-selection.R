test_that("Stage 2 accepts an existing Harmony RNA reduction", {
  skip_if_not_installed("SeuratObject")

  counts <- Matrix::Matrix(
    matrix(
      c(1, 0, 2, 0, 1, 3),
      nrow = 2,
      dimnames = list(c("G1", "G2"), c("C1", "C2", "C3"))
    ),
    sparse = TRUE
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts, assay = "RNA")
  object[["ATAC"]] <- SeuratObject::CreateAssayObject(counts = counts)
  object$condition <- c("A", "A", "B")
  object$cell_type <- c("T", "T", "T")

  harmony <- matrix(
    c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
    nrow = 3,
    dimnames = list(colnames(object), c("harmony_1", "harmony_2"))
  )
  lsi <- matrix(
    c(0.6, 0.5, 0.4, 0.3, 0.2, 0.1),
    nrow = 3,
    dimnames = list(colnames(object), c("LSI_1", "LSI_2"))
  )
  object[["harmony"]] <- SeuratObject::CreateDimReducObject(
    embeddings = harmony,
    key = "harmony_",
    assay = "RNA"
  )
  object[["lsi"]] <- SeuratObject::CreateDimReducObject(
    embeddings = lsi,
    key = "LSI_",
    assay = "ATAC"
  )

  contract <- .rc_condition_metacell_cache_contract(
    object = object,
    condition_col = "condition",
    celltype_col = "cell_type",
    rna_assay = "RNA",
    atac_assay = "ATAC",
    metacell_args = list(
      rna_reduction = "harmony",
      rna_dims = 1:2,
      atac_reduction = "lsi",
      atac_dims = 1:2
    )
  )

  expect_identical(contract$analysis_args$rna_reduction, "harmony")
  expect_identical(contract$analysis_args$rna_dims, 1:2)
  expect_identical(contract$rna_reduction$reduction, "harmony")
  expect_identical(contract$rna_reduction$dims, 1:2)
  expect_true(nzchar(contract$rna_reduction$embedding$values_md5))
})

test_that("Harmony selection is forwarded through the canonical metacell path", {
  low_level <- formals(rc_make_supercell2_metacells)
  expect_true(all(c(
    "rna_reduction", "rna_dims", "atac_reduction", "atac_dims"
  ) %in% names(low_level)))
  expect_identical(low_level$rna_reduction, "pca")
  expect_identical(low_level$atac_reduction, "lsi")

  forwarding <- paste(
    deparse(body(.rc_make_condition_pooled_metacells)),
    collapse = "\n"
  )
  expect_match(
    forwarding,
    "c\\(defaults,\\s*metacell_args\\)",
    perl = TRUE
  )
  expect_false(grepl(
    'c\\([^)]*"rna_reduction"',
    forwarding,
    perl = TRUE
  ))
})

test_that("a missing requested Harmony reduction fails before checkpoint reuse", {
  skip_if_not_installed("SeuratObject")

  counts <- Matrix::Matrix(
    matrix(
      c(1, 0, 2, 0, 1, 3),
      nrow = 2,
      dimnames = list(c("G1", "G2"), c("C1", "C2", "C3"))
    ),
    sparse = TRUE
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts, assay = "RNA")
  object[["ATAC"]] <- SeuratObject::CreateAssayObject(counts = counts)
  object$condition <- c("A", "A", "B")
  object$cell_type <- c("T", "T", "T")
  object[["lsi"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      seq_len(6),
      nrow = 3,
      dimnames = list(colnames(object), c("LSI_1", "LSI_2"))
    ),
    key = "LSI_",
    assay = "ATAC"
  )

  expect_error(
    .rc_condition_metacell_cache_contract(
      object = object,
      condition_col = "condition",
      celltype_col = "cell_type",
      rna_assay = "RNA",
      atac_assay = "ATAC",
      metacell_args = list(
        rna_reduction = "harmony",
        rna_dims = 1:2,
        atac_reduction = "lsi",
        atac_dims = 1:2
      )
    ),
    "Required metacell reduction is absent: `harmony`"
  )
})
