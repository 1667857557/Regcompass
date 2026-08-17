test_that("Stage 2 accepts an existing Harmony RNA reduction", {
  skip_if_not_installed("SeuratObject")
  counts <- Matrix::Matrix(
    matrix(c(1, 0, 2, 0, 1, 3), nrow = 2,
      dimnames = list(c("G1", "G2"), c("C1", "C2", "C3"))),
    sparse = TRUE
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts, assay = "RNA")
  object[["ATAC"]] <- SeuratObject::CreateAssayObject(counts = counts)
  object$condition <- c("A", "A", "B")
  object$cell_type <- c("T", "T", "T")
  harmony <- matrix(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6), nrow = 3,
    dimnames = list(colnames(object), c("harmony_1", "harmony_2")))
  lsi <- matrix(c(0.6, 0.5, 0.4, 0.3, 0.2, 0.1), nrow = 3,
    dimnames = list(colnames(object), c("LSI_1", "LSI_2")))
  object[["harmony"]] <- SeuratObject::CreateDimReducObject(
    embeddings = harmony, key = "harmony_", assay = "RNA"
  )
  object[["lsi"]] <- SeuratObject::CreateDimReducObject(
    embeddings = lsi, key = "LSI_", assay = "ATAC"
  )
  contract <- .rc_condition_metacell_cache_contract(
    object = object, condition_col = "condition", celltype_col = "cell_type",
    rna_assay = "RNA", atac_assay = "ATAC",
    metacell_args = list(
      rna_reduction = "harmony", rna_dims = 1:2,
      atac_reduction = "lsi", atac_dims = 1:2
    )
  )
  expect_identical(contract$analysis_args$rna_reduction, "harmony")
  expect_identical(contract$analysis_args$rna_dims, 1:2)
  expect_identical(contract$analysis_args$gamma, 30L)
  expect_identical(contract$rna_reduction$reduction, "harmony")
  expect_identical(contract$rna_reduction$dims, 1:2)
  expect_true(nzchar(contract$rna_reduction$embedding$values_md5))
  expect_identical(
    contract$native_supercell_api,
    "SCimplify_by_graph_group"
  )
  expect_identical(
    contract$schema_version,
    "regcompass_shared_walktrap_condition_local_repair_cache_v2"
  )
  expect_identical(
    contract$partition_schema_version,
    "shared_walktrap_condition_local_repair_v2"
  )
  expect_identical(
    contract$graph_scope,
    "one_independent_WNN_graph_per_cell_type"
  )
  expect_identical(
    contract$condition_scope,
    "shared_WNN_and_Walktrap_with_condition_gamma_cut_and_local_tree_repair"
  )
  expect_identical(
    contract$membership_split_timing,
    "condition_gamma_cut_then_local_same_condition_hierarchy_repair"
  )
})

test_that("Harmony selection reaches the grouped WNN SuperCell path", {
  builder <- paste(
    deparse(body(.rc_build_grouped_wnn_membership)), collapse = "\n"
  )
  wrapper <- paste(
    deparse(body(.rc_make_condition_celltype_metacells)), collapse = "\n"
  )
  expect_match(builder, "SCimplify_by_graph_group", fixed = TRUE)
  expect_match(
    builder,
    "reduction = list(rna_reduction, atac_reduction)",
    fixed = TRUE
  )
  expect_match(
    builder,
    "dims = list(as.integer(rna_dims), as.integer(atac_dims))",
    fixed = TRUE
  )
  expect_match(wrapper, "rna_reduction = args$rna_reduction", fixed = TRUE)
  expect_match(wrapper, "atac_reduction = args$atac_reduction", fixed = TRUE)
  expect_match(
    wrapper,
    "condition_gamma_cut_then_local_same_condition_hierarchy_repair",
    fixed = TRUE
  )
  expect_false(exists(".rc_native_supercell_membership", inherits = TRUE))
  expect_false(exists(".rc_scale_embedding_block_by_group", inherits = TRUE))
})

test_that("a missing requested Harmony reduction fails before checkpoint reuse", {
  skip_if_not_installed("SeuratObject")
  counts <- Matrix::Matrix(
    matrix(c(1, 0, 2, 0, 1, 3), nrow = 2,
      dimnames = list(c("G1", "G2"), c("C1", "C2", "C3"))),
    sparse = TRUE
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts, assay = "RNA")
  object[["ATAC"]] <- SeuratObject::CreateAssayObject(counts = counts)
  object$condition <- c("A", "A", "B")
  object$cell_type <- c("T", "T", "T")
  object[["lsi"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(seq_len(6), nrow = 3,
      dimnames = list(colnames(object), c("LSI_1", "LSI_2"))),
    key = "LSI_", assay = "ATAC"
  )
  expect_error(
    .rc_condition_metacell_cache_contract(
      object = object, condition_col = "condition", celltype_col = "cell_type",
      rna_assay = "RNA", atac_assay = "ATAC",
      metacell_args = list(
        rna_reduction = "harmony", rna_dims = 1:2,
        atac_reduction = "lsi", atac_dims = 1:2
      )
    ),
    "Required metacell reduction is absent: `harmony`"
  )
})
