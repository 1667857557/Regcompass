test_that("v1.8.0 public workflow is GRN first", {
  text <- paste(deparse(body(rc_run_regcompass)), collapse = "\n")
  stages <- c("rc_regcompass_step_grn", "rc_regcompass_step_metacells", "rc_regcompass_step_meta_modules")
  positions <- vapply(stages, function(x) regexpr(x, text, fixed = TRUE)[[1L]], integer(1))
  expect_true(all(positions > 0L))
  expect_true(positions[[1L]] < positions[[2L]] && positions[[2L]] < positions[[3L]])
  expect_false("inference_unit" %in% names(formals(rc_run_regcompass)))
  expect_identical(eval(formals(rc_run_regcompass)$fragment_files), FALSE)
})

test_that("Pando grouping is condition by cell type", {
  text <- paste(deparse(body(.rc_run_condition_single_cell_grns)), collapse = "\n")
  expect_match(text, "group_cols <- c(condition_col, celltype_col)", fixed = TRUE)
  expect_match(text, "Every condition-by-cell-type Pando GRN must complete successfully", fixed = TRUE)
})

test_that("global union contains biological and local support reactions", {
  artifact <- list(group_id = "condition_pooled", grn_meta_modules = list(
    sample_status = data.frame(status = "ok"), tf_peak_gene_all = data.frame(),
    tf_peak_gene_significant = data.frame(), metabolic_gene_nodes = data.frame(),
    metabolic_gene_edges = data.frame(),
    core_gene_reaction = data.frame(sample_id = "A|T", module_id = "M1", reaction_id = "R1", is_core = TRUE),
    reaction_membership = data.frame(sample_id = "A|T", module_id = "M1", reaction_id = c("R1", "R2")),
    local_completed_reaction_membership = data.frame(sample_id = "A|T", module_id = "M1", reaction_id = c("R1", "R2", "R3")),
    meta_module_summary = data.frame()))
  out <- .rc_merge_stratum_meta_modules(list(artifact))
  expect_setequal(out$global_reaction_membership$reaction_id, c("R1", "R2", "R3"))
  expect_equal(out$global_reaction_membership$inclusion_stage[
    out$global_reaction_membership$reaction_id == "R3"],
    "global_union_local_fastcore_support")
})

test_that("condition metacells have no sample balancing path", {
  text <- paste(deparse(body(.rc_make_condition_pooled_metacells)), collapse = "\n")
  expect_match(text, "condition_x_celltype", fixed = TRUE)
  expect_match(text, 'sample_weighting <- "none"', fixed = TRUE)
  expect_match(text, "gamma <- 75L", fixed = TRUE)
  expect_match(text, "Sample balancing is not part", fixed = TRUE)
})

test_that("condition metacells reject fragment pooling without maps", {
  text <- paste(deparse(body(.rc_make_condition_pooled_metacells)), collapse = "\n")
  expect_match(text, "requires `fragment_files = FALSE`", fixed = TRUE)
})
