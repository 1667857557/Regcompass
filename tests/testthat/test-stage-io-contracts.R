test_that("GRN and metacell groups require bidirectional coverage", {
  grn <- list(group_status = data.frame(
    condition = c("A", "B"), cell_type = c("T", "T"),
    status = c("ok", "ok"), n_cells = c(100L, 120L),
    n_active_edges = c(10L, 0L), stringsAsFactors = FALSE
  ))
  metacells <- data.frame(
    metacell_id = c("A1", "A2", "B1"),
    condition = c("A", "A", "B"), cell_type = "T",
    dominant_celltype_fraction = c(1, 1, 1),
    mixed_celltype_metacell = FALSE,
    stringsAsFactors = FALSE
  )
  coverage <- .rc_validate_grn_metacell_group_coverage(
    grn, metacells, "condition", "cell_type"
  )
  expect_true(all(coverage$coverage_complete))
  expect_true(coverage$has_active_grn_evidence[coverage$condition == "A"])
  expect_false(coverage$has_active_grn_evidence[coverage$condition == "B"])
  expect_equal(coverage$n_metacells[coverage$condition == "A"], 2L)
  expect_error(
    .rc_validate_grn_metacell_group_coverage(
      grn, metacells[metacells$condition == "A", , drop = FALSE],
      "condition", "cell_type"
    ),
    "do not align"
  )
})

test_that("cell-type audit rejects tied memberships", {
  skip_if_not_installed("SeuratObject")
  counts <- Matrix::Matrix(
    matrix(seq_len(8), nrow = 2,
           dimnames = list(c("g1", "g2"), paste0("c", 1:4))),
    sparse = TRUE
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object$cell_type <- c("T", "B", "T", "B")
  pooled <- list(
    membership = data.frame(
      cell_id = paste0("c", 1:4),
      metacell_id = c("mc1", "mc1", "mc2", "mc2"),
      stringsAsFactors = FALSE
    ),
    metacell_meta = data.frame(
      metacell_id = c("mc1", "mc2"), stringsAsFactors = FALSE
    )
  )
  expect_error(
    .rc_assign_metacell_dominant_celltype(pooled, object, "cell_type"),
    "tied cell-type labels"
  )
})

test_that("metacell stage persists required condition-only artifacts", {
  text <- paste(deparse(body(rc_regcompass_step_metacells)), collapse = "\n")
  required <- c(
    "metacell_metadata.tsv.gz", "metacell_membership.tsv.gz",
    "metacell_celltype_composition.tsv.gz",
    "metacell_celltype_summary.tsv.gz", "merged_metacell_object.rds",
    "step_metacells.rds"
  )
  expect_true(all(vapply(required, grepl, logical(1), x = text, fixed = TRUE)))
  expect_false("sample_col" %in% names(formals(rc_regcompass_step_metacells)))
})

test_that("Stage 1 persists shared candidates and bootstrap condition edges", {
  text <- paste(deparse(body(.rc_run_celltype_multitask_grns)), collapse = "\n")
  required <- c(
    "pando_tf_peak_gene_candidates.tsv.gz",
    "pando_tf_peak_gene_global.tsv.gz",
    "pando_tf_peak_gene_condition_all.tsv.gz",
    "pando_tf_peak_gene_significant.tsv.gz",
    "condition_target_genes.tsv.gz", "target_model_diagnostics.tsv.gz",
    "bootstrap_stability_diagnostics.tsv.gz", "pando_celltype_status.tsv.gz",
    "pando_group_status.tsv.gz"
  )
  expect_true(all(vapply(required, grepl, logical(1), x = text, fixed = TRUE)))
  expect_match(text, "regcompass_multitask_grn_v2", fixed = TRUE)
  expect_match(text, "multitask_shared_backbone", fixed = TRUE)
  expect_match(text, "condition_stratified_full_size_nonparametric", fixed = TRUE)
  expect_match(text, "pando_grn_design_v2", fixed = TRUE)
})

test_that("Stage 3 persists bootstrap-supported genes and complete cores", {
  text <- paste(deparse(body(.rc_build_condition_meta_modules)), collapse = "\n")
  required <- c(
    "supported_metabolic_genes.tsv.gz", "core_gene_reaction.tsv.gz",
    "meta_module_reactions.tsv.gz", "condition_meta_modules.rds"
  )
  expect_true(all(vapply(required, grepl, logical(1), x = text, fixed = TRUE)))
  expect_match(text, ".rc_summarize_supported_metabolic_genes", fixed = TRUE)
  expect_match(text, "group_id", fixed = TRUE)
  expect_false(grepl("rc_project_metabolic_grn", text, fixed = TRUE))
})

test_that("Layer 1 records exact SuperCell2 label provenance", {
  body_text <- paste(
    deparse(body(.rc_build_condition_pooled_layer1)), collapse = "\n"
  )
  step_text <- paste(deparse(body(rc_regcompass_step_layer1)), collapse = "\n")
  expect_match(
    body_text, "regcompass_condition_only_layer1_v5_supercell_label",
    fixed = TRUE
  )
  expect_match(
    body_text, "condition_stratified_supercell2_label_exact_celltype",
    fixed = TRUE
  )
  expect_match(body_text, "Layer 1 requires SuperCell2 label-pure", fixed = TRUE)
  expect_match(body_text, "and_method = gpr_and_method", fixed = TRUE)
  expect_false("sample_col" %in% names(formals(.rc_build_condition_pooled_layer1)))
  expect_match(step_text, "regcompass_layer1_step", fixed = TRUE)
})

test_that("Layer 2 and compact final results validate provenance", {
  layer2_text <- paste(deparse(body(rc_regcompass_step_layer2)), collapse = "\n")
  result_text <- paste(deparse(body(rc_regcompass_step_results)), collapse = "\n")
  expect_match(layer2_text, ".rc_validate_layer1_stage", fixed = TRUE)
  expect_match(layer2_text, "regcompass_layer2_step", fixed = TRUE)
  expect_match(result_text, ".rc_validate_layer2_stage", fixed = TRUE)
  expect_match(
    result_text, "regcompass_compact_multitask_result_v3", fixed = TRUE
  )
  expect_match(result_text, 'version = "1.8.9"', fixed = TRUE)
  expect_match(result_text, "active_regulatory_edges", fixed = TRUE)
  expect_match(result_text, "condition_target_genes", fixed = TRUE)
  expect_match(result_text, "core_reactions", fixed = TRUE)
  expect_match(result_text, "table_manifest", fixed = TRUE)
  expect_match(result_text, "stage_provenance", fixed = TRUE)
  expect_false(grepl("condition_grn_meta_modules", result_text, fixed = TRUE))
  expect_false(grepl("merged_grn_meta_modules", result_text, fixed = TRUE))
})

test_that("stage validators reject reordered or mismatched units", {
  params <- list(
    condition_col = "condition", celltype_col = "cell_type",
    rna_assay = "RNA", atac_assay = "ATAC"
  )
  layer1 <- list(
    reaction_expression = matrix(
      1, nrow = 1, ncol = 2,
      dimnames = list("R1", c("U1", "U2"))
    ),
    unit_meta = data.frame(
      pool_id = c("U2", "U1"), stringsAsFactors = FALSE
    ),
    workflow_params = params,
    gem_fingerprint = "x"
  )
  class(layer1) <- c("regcompass_layer1_step", "list")
  expect_error(.rc_validate_layer1_stage(layer1), "not identically ordered")
})
