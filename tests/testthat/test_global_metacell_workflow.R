test_that("v1.7.0 public workflow is condition pooled", {
  wrapper_text <- paste(deparse(body(rc_run_regcompass)), collapse = "\n")
  base_text <- paste(deparse(body(.rc_run_regcompass_uncorrected_metadata)),
                     collapse = "\n")
  expect_match(base_text, ".rc_make_condition_pooled_metacells", fixed = TRUE)
  expect_match(base_text, ".rc_run_condition_pando_modules", fixed = TRUE)
  expect_match(base_text, "missing_conditions", fixed = TRUE)
  pando_text <- paste(deparse(body(.rc_run_condition_pando_modules)), collapse = "\n")
  expect_match(pando_text, "Every condition-by-cell-type Pando GRN must complete successfully", fixed = TRUE)
  expect_match(base_text, ".rc_merge_stratum_meta_modules", fixed = TRUE)
  expect_match(base_text, ".rc_build_condition_pooled_layer1", fixed = TRUE)
  expect_match(base_text, "04_model_cache", fixed = TRUE)
  expect_match(wrapper_text, ".rc_condition_pool_design_summary", fixed = TRUE)
  expect_false("inference_unit" %in% names(formals(rc_run_regcompass)))
  expect_identical(eval(formals(rc_run_regcompass)$fragment_files), FALSE)
})

test_that("condition grouping excludes biological sample", {
  expect_identical(
    .rc_condition_group_cols("condition", "cell_type"),
    c("condition", "cell_type")
  )
  expect_error(
    .rc_condition_group_cols("condition", "condition"),
    "must be distinct"
  )
})

test_that("global union contains biological and local support reactions", {
  artifact <- list(
    group_id = "condition_pooled",
    grn_meta_modules = list(
      sample_status = data.frame(status = "ok"),
      tf_peak_gene_all = data.frame(),
      tf_peak_gene_significant = data.frame(),
      metabolic_gene_nodes = data.frame(),
      metabolic_gene_edges = data.frame(),
      core_gene_reaction = data.frame(
        sample_id = "A|T",
        module_id = "M1",
        reaction_id = "R1",
        is_core = TRUE
      ),
      reaction_membership = data.frame(
        sample_id = "A|T",
        module_id = "M1",
        reaction_id = c("R1", "R2")
      ),
      local_completed_reaction_membership = data.frame(
        sample_id = "A|T",
        module_id = "M1",
        reaction_id = c("R1", "R2", "R3")
      ),
      meta_module_summary = data.frame()
    )
  )
  out <- .rc_merge_stratum_meta_modules(list(artifact))
  expect_setequal(out$global_reaction_membership$reaction_id, c("R1", "R2", "R3"))
  expect_equal(
    out$global_reaction_membership$inclusion_stage[
      out$global_reaction_membership$reaction_id == "R3"
    ],
    "global_union_local_fastcore_support"
  )
})

test_that("condition-pooled workflow balances samples before base pooling", {
  wrapper_text <- paste(
    deparse(body(.rc_make_condition_pooled_metacells)),
    collapse = "\n"
  )
  expect_match(wrapper_text, ".rc_balance_condition_celltype_cells", fixed = TRUE)
  expect_match(wrapper_text, "sample_balance_diagnostics.tsv.gz", fixed = TRUE)
  expect_match(wrapper_text, "sample_balance %||% TRUE", fixed = TRUE)
})

test_that("condition-pooled workflow rejects fragment pooling without maps", {
  text <- paste(
    deparse(body(.rc_make_condition_pooled_metacells_unbalanced)),
    collapse = "\n"
  )
  expect_match(text, "fragment_files = FALSE", fixed = TRUE)
  expect_match(text, "explicit per-file barcode map", fixed = TRUE)
})
