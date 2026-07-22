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

test_that("metacell construction is condition-only without sample balancing", {
  text <- paste(deparse(body(.rc_make_condition_pooled_metacells)), collapse = "\n")
  expect_match(text, 'object@meta.data[[internal_celltype_col]] <- "all_celltypes"', fixed = TRUE)
  expect_match(text, 'pooling_scope <- "condition_only"', fixed = TRUE)
  expect_match(text, 'sample_weighting <- "none"', fixed = TRUE)
  expect_match(text, "metacell_grouping = condition_col", fixed = TRUE)
  expect_match(text, "gamma <- 75L", fixed = TRUE)
  expect_match(text, "Sample balancing is not part", fixed = TRUE)
  expect_match(text, "label_col = label_col", fixed = TRUE)
})

test_that("public metacell workflows expose a pre-aggregation annotation label", {
  step_formals <- formals(rc_regcompass_step_metacells)
  run_formals <- formals(rc_run_regcompass)

  expect_identical(deparse(step_formals$label_col), "celltype_col")
  expect_identical(deparse(run_formals$metacell_label_col), "celltype_col")

  step_text <- paste(deparse(body(rc_regcompass_step_metacells)), collapse = "\n")
  run_text <- paste(deparse(body(rc_run_regcompass)), collapse = "\n")
  expect_match(step_text, "label_col = label_col", fixed = TRUE)
  expect_match(run_text, "label_col = metacell_label_col", fixed = TRUE)
})

test_that("dominant cell type is assigned after condition-only metacells", {
  skip_if_not_installed("SeuratObject")
  counts <- Matrix::Matrix(matrix(1, nrow = 1, ncol = 6,
    dimnames = list("g1", paste0("c", 1:6))), sparse = TRUE)
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object$cell_type <- c("T", "T", "B", "B", "B", "T")
  pooled <- list(
    membership = data.frame(
      cell_id = paste0("c", 1:6),
      metacell_id = c("m1", "m1", "m1", "m2", "m2", "m2")
    ),
    metacell_meta = data.frame(metacell_id = c("m1", "m2"))
  )
  out <- .rc_assign_metacell_dominant_celltype(pooled, object, "cell_type")
  expect_identical(out$metacell_meta$cell_type, c("T", "B"))
  expect_equal(out$metacell_meta$dominant_celltype_fraction, c(2/3, 2/3))
  expect_true(all(out$metacell_meta$mixed_celltype_metacell))
  expect_false(any(out$metacell_meta$dominant_celltype_tied))
})

test_that("condition metacells reject fragment pooling without maps", {
  text <- paste(deparse(body(.rc_make_condition_pooled_metacells)), collapse = "\n")
  expect_match(text, "requires `fragment_files = FALSE`", fixed = TRUE)
})
