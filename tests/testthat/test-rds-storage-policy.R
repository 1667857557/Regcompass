test_that("duplicate stage sidecar RDS files are not written", {
  writer <- getFromNamespace("saveRDS", "RegCompassR")
  outdir <- tempfile("regcompass-rds-policy-")
  dir.create(outdir, recursive = TRUE)

  duplicate_files <- c(
    "single_cell_grn.rds",
    "pando_condition_grn_fits.rds",
    "condition_grn_fit.rds",
    "standard_pando_T_cell.rds",
    "rna_counts.rds",
    "atac_counts.rds",
    "metacell_object.rds",
    "merged_metacell_object.rds",
    "condition_metacell_cache_contract.rds",
    "merged_meta_modules.rds",
    "strict_score_matrix.rds",
    "strict_penalty_matrix.rds",
    "vmax_matrix.rds",
    "feasible_matrix.rds",
    "evaluated_matrix.rds",
    "penalty_components.rds",
    "run_parameters.rds",
    "model_completion_contract.rds",
    "step_comparison.rds",
    "00_model_info.rds",
    "00_medium_scenarios.rds"
  )

  for (name in duplicate_files) {
    path <- file.path(outdir, name)
    writer(list(payload = rep(1, 10)), path)
    expect_false(file.exists(path), info = name)
  }
})

test_that("externalized Stage 3 and structural model caches remain writable", {
  writer <- getFromNamespace("saveRDS", "RegCompassR")
  outdir <- tempfile("regcompass-rds-essential-")
  dir.create(outdir, recursive = TRUE)

  condition_modules <- file.path(outdir, "condition_meta_modules.rds")
  union_model <- file.path(
    outdir,
    "union_gem__celltype_545f63656c6c__medium_62617365.rds"
  )

  writer(list(reaction_membership = data.frame(reaction_id = "R1")),
         condition_modules)
  writer(list(S = matrix(1, 1, 1)), union_model)

  expect_true(file.exists(condition_modules))
  expect_true(file.exists(union_model))
})

test_that("Stage 1 storage removes exact duplicate GRN aliases", {
  writer <- getFromNamespace("saveRDS", "RegCompassR")
  outdir <- tempfile("stage1-storage-")
  dir.create(outdir, recursive = TRUE)
  path <- file.path(outdir, "step_grn.rds")

  canonical_object <- list(marker = "GRNData")
  x <- list(
    grn_result = list(
      pando_grn_data = canonical_object,
      pando_grn_data_by_cell_type = list(T_cell = canonical_object),
      tf_peak_gene_condition_all = data.frame(edge_id = "E1"),
      tf_peak_gene_condition = data.frame(edge_id = "E1"),
      tf_peak_gene_condition_effect_all = data.frame(edge_id = "E1"),
      tf_peak_gene_condition_effect = data.frame(edge_id = "E1")
    )
  )

  writer(x, path)
  stored <- base::readRDS(path)

  expect_null(stored$grn_result$pando_grn_data)
  expect_equal(names(stored$grn_result$pando_grn_data_by_cell_type), "T_cell")
  expect_true(is.data.frame(stored$grn_result$tf_peak_gene_condition_all))
  expect_true(is.data.frame(stored$grn_result$tf_peak_gene_condition))
  expect_null(stored$grn_result$tf_peak_gene_condition_effect_all)
  expect_null(stored$grn_result$tf_peak_gene_condition_effect)
})

test_that("Stage 2 storage keeps one canonical metacell object", {
  writer <- getFromNamespace("saveRDS", "RegCompassR")
  outdir <- tempfile("stage2-storage-")
  dir.create(outdir, recursive = TRUE)
  path <- file.path(outdir, "step_metacells.rds")

  canonical <- list(id = "canonical_metacell_object", assay_payload = 1:20)
  x <- list(
    pooled = list(
      metacell_objects = list(grouped_wnn = canonical),
      metacell_object = canonical,
      rna_counts = matrix(1, 2, 2),
      atac_counts = matrix(2, 2, 2),
      metacell_meta = data.frame(metacell_id = "MC1"),
      membership = data.frame(cell_id = "C1", metacell_id = "MC1"),
      cache_contract = list(schema_version = "contract"),
      input_design = list(graph_method = "multimodal_WNN")
    ),
    metacell_object = canonical
  )

  writer(x, path)
  stored <- base::readRDS(path)

  expect_equal(stored$metacell_object$id, "canonical_metacell_object")
  expect_null(stored$pooled$metacell_objects)
  expect_null(stored$pooled$metacell_object)
  expect_null(stored$pooled$rna_counts)
  expect_null(stored$pooled$atac_counts)
  expect_true(is.data.frame(stored$pooled$metacell_meta))
  expect_true(is.data.frame(stored$pooled$membership))
})

test_that("Stage 5 storage drops the duplicate RNA-only result route", {
  writer <- getFromNamespace("saveRDS", "RegCompassR")
  outdir <- tempfile("stage5-storage-")
  dir.create(outdir, recursive = TRUE)
  path <- file.path(outdir, "step_layer2.rds")

  x <- list(
    penalty = matrix(1, 2, 2),
    penalty_rna_only = matrix(2, 2, 2),
    score_rna_only = matrix(3, 2, 2),
    comparison_paths = list(
      rna_only = list(
        penalty = matrix(2, 2, 2),
        score = matrix(3, 2, 2),
        penalty_components = matrix(4, 20, 20)
      )
    )
  )

  writer(x, path)
  stored <- base::readRDS(path)

  expect_true(is.matrix(stored$penalty_rna_only))
  expect_true(is.matrix(stored$score_rna_only))
  expect_null(stored$comparison_paths)
})

test_that("saved final result omits complete upstream stage payloads", {
  writer <- getFromNamespace("saveRDS", "RegCompassR")
  outdir <- tempfile("final-storage-")
  dir.create(outdir, recursive = TRUE)
  path <- file.path(outdir, "regcompass_result.rds")

  x <- list(
    schema_version = "regcompass_regulatory_metabolic_result_v3",
    grn = list(big = rep(1, 100)),
    metacells = list(big = rep(2, 100)),
    layer1 = list(big = rep(3, 100)),
    condition_grn_meta_modules = list(big = rep(4, 100)),
    merged_grn_meta_modules = list(big = rep(5, 100)),
    grn_meta_modules = list(big = rep(6, 100)),
    microcompass = list(big = rep(7, 100)),
    reaction_catalog = data.frame(reaction_id = "R1"),
    reaction_evidence = data.frame(reaction_id = "R1"),
    reaction_comparison_by_metacell = data.frame(reaction_id = "R1"),
    reaction_ranking = data.frame(reaction_id = "R1"),
    condition_summary = data.frame(reaction_id = "R1"),
    condition_contrast = data.frame(),
    rna_only_control_summary = data.frame(reaction_id = "R1"),
    rna_only_control_contrast = data.frame()
  )

  writer(x, path)
  stored <- base::readRDS(path)

  expect_null(stored$grn)
  expect_null(stored$metacells)
  expect_null(stored$layer1)
  expect_null(stored$condition_grn_meta_modules)
  expect_null(stored$merged_grn_meta_modules)
  expect_null(stored$grn_meta_modules)
  expect_null(stored$microcompass)
  expect_true(is.data.frame(stored$reaction_catalog))
  expect_true(is.data.frame(stored$reaction_comparison_by_metacell))
  expect_identical(
    stored$storage_contract$schema_version,
    "regcompass_compact_final_result_storage_v1"
  )
  expect_false(stored$storage_contract$upstream_stage_payloads_embedded)
})

test_that("one-shot root final result is not duplicated after Stage 6", {
  writer <- getFromNamespace("saveRDS", "RegCompassR")
  root <- tempfile("one-shot-result-")
  dir.create(file.path(root, "06_results"), recursive = TRUE)
  canonical <- file.path(root, "06_results", "regcompass_result.rds")
  duplicate <- file.path(root, "regcompass_result.rds")

  writer(list(reaction_catalog = data.frame(reaction_id = "R1")), canonical)
  writer(list(reaction_catalog = data.frame(reaction_id = "R1")), duplicate)

  expect_true(file.exists(canonical))
  expect_false(file.exists(duplicate))
})
