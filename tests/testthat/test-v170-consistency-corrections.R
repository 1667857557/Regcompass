test_that("cell-type-shared TF-IDF is reused across conditions", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Signac")

  cells <- c("T_A", "T_B", "B_A", "B_B")
  rna <- Matrix::Matrix(
    matrix(
      c(5, 4, 3, 2,
        2, 3, 4, 5),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("G1", "G2"), cells)
    ),
    sparse = TRUE
  )
  atac <- Matrix::Matrix(
    matrix(
      c(5, 1, 2, 2,
        1, 5, 2, 2,
        1, 1, 5, 1),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("p1", "p2", "p3"), cells)
    ),
    sparse = TRUE
  )
  object <- Seurat::CreateSeuratObject(counts = rna, assay = "RNA")
  object[["ATAC"]] <- SeuratObject::CreateAssayObject(counts = atac)
  object$condition <- c("A", "B", "A", "B")
  object$cell_type <- c("T", "T", "B", "B")

  normalized <- .rc_apply_celltype_shared_tfidf(
    object,
    celltype_col = "cell_type",
    atac_assay = "ATAC"
  )
  observed <- .rc_pando_assay_data(normalized, "ATAC")
  expected_t <- Signac::RunTFIDF(
    atac[, c("T_A", "T_B"), drop = FALSE],
    verbose = FALSE
  )
  expected_b <- Signac::RunTFIDF(
    atac[, c("B_A", "B_B"), drop = FALSE],
    verbose = FALSE
  )

  expect_equal(
    as.matrix(observed[, c("T_A", "T_B"), drop = FALSE]),
    as.matrix(expected_t)
  )
  expect_equal(
    as.matrix(observed[, c("B_A", "B_B"), drop = FALSE]),
    as.matrix(expected_b)
  )
  expect_identical(
    normalized@misc$regcompass_atac_normalization$scope,
    "cell_type_across_conditions"
  )
  expect_identical(
    normalized@misc$regcompass_atac_normalization$n_metacells_by_celltype,
    c(B = 2L, T = 2L)
  )
})

test_that("custom cell-type metadata is inferred from pooled design", {
  pooled <- list(
    input_design = list(
      condition_celltype_sample_count = data.frame(
        treatment = "A",
        epithelial_or_stem = "epithelial_like",
        n_biological_samples = 2L,
        check.names = FALSE
      )
    ),
    metacell_meta = data.frame(epithelial_or_stem = "epithelial_like")
  )
  expect_identical(
    .rc_celltype_col_from_pooled(pooled),
    "epithelial_or_stem"
  )
})

test_that("Pando reuses Step 1 normalized data without group renormalization", {
  body_text <- paste(deparse(body(rc_run_pando_meta_modules)), collapse = "\n")
  expect_false(grepl("Signac::RunTFIDF", body_text, fixed = TRUE))
  expect_false(grepl("Seurat::NormalizeData", body_text, fixed = TRUE))
  expect_match(body_text, ".rc_require_normalized_assay", fixed = TRUE)
  expect_match(body_text, "cell_type_across_conditions", fixed = TRUE)
})

test_that("targets without finite Pando R-squared are untrusted", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")

  counts <- Matrix::Matrix(
    matrix(
      c(2, 3, 1, 1),
      nrow = 2,
      dimnames = list(c("G1", "G2"), c("mc1", "mc2"))
    ),
    sparse = TRUE
  )
  object <- Seurat::CreateSeuratObject(counts = counts, assay = "RNA")
  edges <- data.frame(
    target = "GENE1",
    region = "p1",
    tf = "TF1",
    estimate = 1,
    rsq = NA_real_,
    condition = "A",
    cell_type = "T",
    stringsAsFactors = FALSE
  )
  unit_meta <- data.frame(
    pool_id = c("mc1", "mc2"),
    condition = "A",
    cell_type = "T",
    stringsAsFactors = FALSE
  )

  modifier <- .rc_condition_gene_regulatory_modifier(
    significant_edges = edges,
    object = object,
    unit_meta = unit_meta,
    target_genes = "gene1"
  )

  expect_true(all(modifier == 0))
  expect_match(
    attr(modifier, "reliability_policy"),
    "reliability zero",
    fixed = TRUE
  )
})

test_that("meta-module GEM completion provenance describes both stages", {
  meta_module <- .rc_feasibility_completion_metadata("meta_module_gem")
  expect_identical(
    meta_module$feasibility_completion,
    "local_unconstrained_fastcore_then_global_union_medium_specific_fastcore"
  )
  expect_match(meta_module$feasibility_completion_stages$local, "unconstrained")
  expect_match(meta_module$feasibility_completion_stages$global, "medium")

  full_gem <- .rc_feasibility_completion_metadata("full_gem")
  expect_identical(
    full_gem$feasibility_completion,
    "not_applicable_full_gem"
  )
})
