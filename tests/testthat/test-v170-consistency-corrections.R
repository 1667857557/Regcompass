test_that("cell-type-shared TF-IDF is computed across available cells", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Signac")
  cells <- c("T_A", "T_B", "B_A", "B_B")
  rna <- Matrix::Matrix(
    matrix(c(5, 4, 3, 2, 2, 3, 4, 5), nrow = 2, byrow = TRUE,
      dimnames = list(c("G1", "G2"), cells)), sparse = TRUE
  )
  atac <- Matrix::Matrix(
    matrix(c(5, 1, 2, 2, 1, 5, 2, 2, 1, 1, 5, 1), nrow = 3,
      byrow = TRUE, dimnames = list(c("p1", "p2", "p3"), cells)),
    sparse = TRUE
  )
  object <- Seurat::CreateSeuratObject(counts = rna, assay = "RNA")
  object[["ATAC"]] <- SeuratObject::CreateAssayObject(counts = atac)
  object$condition <- c("A", "B", "A", "B")
  object$cell_type <- c("T", "T", "B", "B")
  normalized <- .rc_apply_celltype_shared_tfidf(
    object, "cell_type", "ATAC"
  )
  observed <- .rc_pando_assay_data(normalized, "ATAC")
  expect_equal(
    as.matrix(observed[, c("T_A", "T_B"), drop = FALSE]),
    as.matrix(Signac::RunTFIDF(
      atac[, c("T_A", "T_B"), drop = FALSE], verbose = FALSE
    ))
  )
  expect_equal(
    as.matrix(observed[, c("B_A", "B_B"), drop = FALSE]),
    as.matrix(Signac::RunTFIDF(
      atac[, c("B_A", "B_B"), drop = FALSE], verbose = FALSE
    ))
  )
  expect_identical(
    normalized@misc$regcompass_atac_normalization$scope,
    "cell_type_all_available_cells"
  )
  expect_setequal(
    normalized@misc$regcompass_atac_normalization$n_units_by_celltype,
    c(B = 2L, T = 2L)
  )
})

test_that("single-cell Pando reuses shared normalized data", {
  condition_impl <- paste(
    deparse(body(.rc_fit_condition_grns_by_cell_type)), collapse = "\n"
  )
  standard_impl <- paste(
    deparse(body(.rc_fit_standard_pando_by_cell_type)), collapse = "\n"
  )
  expect_false(grepl("Signac::RunTFIDF", condition_impl, fixed = TRUE))
  expect_false(grepl("Seurat::NormalizeData", condition_impl, fixed = TRUE))
  expect_false(grepl("Signac::RunTFIDF", standard_impl, fixed = TRUE))
  expect_match(
    paste(deparse(formals(.rc_fit_condition_grns_by_cell_type)$pando_infer_args), collapse = " "),
    "peak_cor = 0", fixed = TRUE
  )
})

test_that("within-cell-type OOF controls condition Layer 1 reliability", {
  projection_text <- paste(
    deparse(body(.rc_condition_pando_projection)), collapse = "\n"
  )
  layer1_text <- paste(
    deparse(body(.rc_cell_first_projection_layer1)), collapse = "\n"
  )
  expect_false(exists(
    ".rc_condition_gene_regulatory_modifier", inherits = TRUE
  ))
  expect_match(projection_text, "q[!available] <- NA_real_", fixed = TRUE)
  expect_match(layer1_text, "gene_regulatory_reliability_available", fixed = TRUE)
  expect_match(
    projection_text, "outer_condition_stratified_cell_oof", fixed = TRUE
  )
})

test_that("feasibility completion is global and medium specific", {
  stage3 <- paste(deparse(body(rc_regcompass_step_meta_modules)), collapse = "\n")
  stage5 <- paste(deparse(body(rc_regcompass_step_layer2)), collapse = "\n")
  cache <- paste(
    deparse(body(.rc_build_medium_specific_union_gem_cache)), collapse = "\n"
  )
  expect_false(grepl(".rc_complete_medium_union_gem", stage3, fixed = TRUE))
  expect_match(stage5, "one union GEM", fixed = TRUE)
  expect_match(cache, ".rc_complete_medium_union_gem", fixed = TRUE)
  expect_match(
    cache, "single_global_fastcore_after_meta_module_merge", fixed = TRUE
  )
})
