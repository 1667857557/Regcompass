test_that("GRN and metacell groups require bidirectional coverage", {
  grn <- list(sample_status = data.frame(
    condition = c("A", "B"),
    cell_type = c("T", "T"),
    status = c("ok", "ok"),
    n_cells = c(100L, 120L),
    n_significant_edges = c(10L, 0L),
    stringsAsFactors = FALSE
  ))
  metacells <- data.frame(
    metacell_id = c("A1", "A2", "B1"),
    condition = c("A", "A", "B"),
    cell_type = "T",
    dominant_celltype_fraction = c(1, 0.8, 1),
    mixed_celltype_metacell = c(FALSE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  coverage <- .rc_validate_grn_metacell_group_coverage(
    grn, metacells, "condition", "cell_type"
  )
  expect_true(all(coverage$coverage_complete))
  expect_true(coverage$has_significant_pando_evidence[coverage$condition == "A"])
  expect_false(coverage$has_significant_pando_evidence[coverage$condition == "B"])
  expect_equal(coverage$n_metacells[coverage$condition == "A"], 2L)
  expect_equal(
    coverage$n_mixed_celltype_metacells[coverage$condition == "A"], 1L
  )
  expect_error(
    .rc_validate_grn_metacell_group_coverage(
      grn,
      metacells[metacells$condition == "A", , drop = FALSE],
      "condition", "cell_type"
    ),
    "do not align"
  )
})

test_that("condition-only metacells reject tied dominant cell types", {
  skip_if_not_installed("SeuratObject")
  counts <- Matrix::Matrix(
    matrix(
      seq_len(8), nrow = 2,
      dimnames = list(c("g1", "g2"), paste0("c", 1:4))
    ),
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
      metacell_id = c("mc1", "mc2"),
      stringsAsFactors = FALSE
    )
  )
  expect_error(
    .rc_assign_metacell_dominant_celltype(pooled, object, "cell_type"),
    "tied dominant cell types"
  )
})

test_that("metacell stage persists required artifacts", {
  text <- paste(deparse(body(rc_regcompass_step_metacells)), collapse = "\n")
  required <- c(
    "metacell_metadata.tsv.gz",
    "metacell_membership.tsv.gz",
    "metacell_celltype_composition.tsv.gz",
    "metacell_celltype_summary.tsv.gz",
    "merged_metacell_object.rds",
    "step_metacells.rds"
  )
  expect_true(all(vapply(required, grepl, logical(1), x = text, fixed = TRUE)))
})

test_that("Stage 3 persists supported genes and core reactions", {
  text <- paste(
    deparse(body(.rc_build_condition_meta_modules)), collapse = "\n"
  )
  required <- c(
    "supported_metabolic_genes.tsv.gz",
    "core_gene_reaction.tsv.gz",
    "meta_module_reactions.tsv.gz",
    "condition_meta_modules.rds"
  )
  expect_true(all(vapply(required, grepl, logical(1), x = text, fixed = TRUE)))
  expect_match(text, ".rc_summarize_supported_metabolic_genes", fixed = TRUE)
  expect_false(grepl("rc_project_metabolic_grn", text, fixed = TRUE))
  expect_false("expansion_mode" %in% names(formals(rc_expand_meta_module_reactions)))
  expect_false("max_iterations" %in% names(formals(rc_expand_meta_module_reactions)))
})

test_that("Layer 1 uses the canonical schema and stage class", {
  body_text <- paste(
    deparse(body(.rc_build_condition_pooled_layer1)), collapse = "\n"
  )
  step_text <- paste(deparse(body(rc_regcompass_step_layer1)), collapse = "\n")
  expect_match(
    body_text,
    "regcompass_condition_only_layer1_v3",
    fixed = TRUE
  )
  expect_match(
    body_text,
    "condition_only_metacell_with_posthoc_celltype",
    fixed = TRUE
  )
  expect_match(body_text, "and_method = gpr_and_method", fixed = TRUE)
  expect_match(step_text, "regcompass_layer1_step", fixed = TRUE)
  expect_match(step_text, "gem_fingerprint", fixed = TRUE)
  expect_match(step_text, "workflow_params", fixed = TRUE)
  expect_identical(
    eval(formals(rc_regcompass_step_layer1)$gpr_and_method),
    c("min", "median", "mean")
  )
})

test_that("Layer 2 and final results validate upstream provenance", {
  layer2_text <- paste(deparse(body(rc_regcompass_step_layer2)), collapse = "\n")
  result_text <- paste(deparse(body(rc_regcompass_step_results)), collapse = "\n")
  expect_match(layer2_text, ".rc_validate_layer1_stage", fixed = TRUE)
  expect_match(layer2_text, "regcompass_layer2_step", fixed = TRUE)
  expect_match(layer2_text, "source_core_reactions", fixed = TRUE)
  expect_match(result_text, ".rc_validate_layer2_stage", fixed = TRUE)
  expect_match(
    result_text,
    "regcompass_significant_pando_targets_v1",
    fixed = TRUE
  )
  expect_match(result_text, 'version = "1.8.4"', fixed = TRUE)
  expect_match(result_text, "condition_grn_meta_modules", fixed = TRUE)
  expect_match(result_text, "merged_grn_meta_modules", fixed = TRUE)
  expect_match(result_text, "supported_metabolic_genes", fixed = TRUE)
  expect_match(result_text, "reaction_catalog", fixed = TRUE)
  expect_match(result_text, "reaction_evidence", fixed = TRUE)
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
      pool_id = c("U2", "U1"),
      stringsAsFactors = FALSE
    ),
    workflow_params = params,
    gem_fingerprint = "x"
  )
  class(layer1) <- c("regcompass_layer1_step", "list")
  expect_error(.rc_validate_layer1_stage(layer1), "not identically ordered")
})
