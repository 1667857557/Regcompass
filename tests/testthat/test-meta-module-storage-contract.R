test_that("condition meta-modules retain only Stage 3-derived payloads", {
  S <- diag(2)
  dimnames(S) <- list(c("M1", "M2"), c("R1", "R2"))
  gem <- rc_make_gem(
    S,
    lb = c(0, 0),
    ub = c(1000, 1000),
    reaction_meta = data.frame(
      reaction_id = c("R1", "R2"),
      subsystem = c("A", "A"),
      metabolic_module = c("A", "A"),
      role = "internal",
      role_source = "test",
      stringsAsFactors = FALSE
    )
  )
  gem$gpr_table <- data.frame(
    reaction_id = "R1",
    and_group_id = "1",
    gene = "G1",
    stringsAsFactors = FALSE
  )
  grn_result <- list(
    group_cols = c("condition", "cell_type"),
    tf_peak_gene_condition = data.frame(
      group_id = "C1|T",
      condition = "C1",
      cell_type = "T",
      tf = "TF1",
      target = "G1",
      region = "chr1-1-2",
      estimate = 1,
      padj = 0.01,
      rsq = 0.5,
      stringsAsFactors = FALSE
    ),
    tf_peak_gene_condition_all = data.frame(stage1 = 1),
    condition_grn_fits = list(stage1 = "do-not-copy"),
    pando_grn_data = list(stage1 = "do-not-copy")
  )
  outdir <- tempfile("regcompass-meta-modules-")
  dir.create(outdir)

  modules <- .rc_build_condition_meta_modules(grn_result, gem, outdir)

  expect_identical(
    modules$schema_version,
    "regcompass_condition_meta_modules_v2"
  )
  expect_true(all(c(
    "supported_metabolic_genes", "core_gene_reaction",
    "reaction_membership", "meta_module_summary", "crossref_maps",
    "core_definition", "expansion_definition", "analysis_group_unit",
    "feasibility_completion"
  ) %in% names(modules)))
  expect_false(any(c(
    "pando_grn_data", "condition_grn_fits",
    "tf_peak_gene_condition", "tf_peak_gene_condition_all",
    "tf_peak_gene_condition_effect", "tf_peak_gene_condition_effect_all",
    "biological_reaction_membership"
  ) %in% names(modules)))

  merged <- .rc_merge_meta_modules_by_cell_type(
    modules,
    celltype_col = "cell_type",
    condition_col = "condition"
  )
  expect_identical(merged$cell_type_catalogues$T$source_conditions, "C1")
  expect_identical(merged$cell_type_catalogues$T$source_group_ids, "C1|T")
  expect_false(any(c(
    "condition_fit_status", "tf_peak_gene_condition_all",
    "tf_peak_gene_condition"
  ) %in% names(merged)))
})
