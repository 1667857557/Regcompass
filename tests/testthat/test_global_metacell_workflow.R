test_that("v1.8.8 public workflow is GRN first and condition only", {
  text <- paste(deparse(body(rc_run_regcompass)), collapse = "\n")
  stages <- c(
    "rc_regcompass_step_grn",
    "rc_regcompass_step_metacells",
    "rc_regcompass_step_meta_modules"
  )
  positions <- vapply(
    stages,
    function(x) regexpr(x, text, fixed = TRUE)[[1L]],
    integer(1)
  )
  expect_true(all(positions > 0L))
  expect_true(
    positions[[1L]] < positions[[2L]] &&
      positions[[2L]] < positions[[3L]]
  )
  expect_false("sample_col" %in% names(formals(rc_run_regcompass)))
  expect_false("sample_col" %in% names(formals(rc_regcompass_step_grn)))
  expect_false("sample_col" %in% names(formals(rc_regcompass_step_metacells)))
  expect_identical(eval(formals(rc_run_regcompass)$fragment_files), FALSE)
  expect_true("meta_module_args" %in% names(formals(rc_run_regcompass)))
  expect_true("multitask_args" %in% names(formals(rc_run_regcompass)))
  expect_identical(
    eval(formals(rc_run_regcompass)$grn_mode),
    c("multitask_shared_backbone", "legacy_condition_pando")
  )
})

test_that("canonical Pando background is shared by cell type", {
  text <- paste(
    deparse(body(.rc_run_celltype_multitask_grns)),
    collapse = "\n"
  )
  expect_match(text, "celltypes <- sort(unique", fixed = TRUE)
  expect_match(text, "Pando::prepare_grn_design", fixed = TRUE)
  expect_match(text, ".rc_fit_multitask_celltype_grn", fixed = TRUE)
  expect_match(text, "edge_universe_id", fixed = TRUE)
  expect_match(text, "pando_grn_design_v2", fixed = TRUE)
  expect_match(
    text,
    "Every cell-type multitask GRN must complete successfully",
    fixed = TRUE
  )
})

test_that("legacy Pando grouping remains available explicitly", {
  text <- paste(
    deparse(body(.rc_run_condition_single_cell_grns)),
    collapse = "\n"
  )
  expect_match(
    text,
    "group_cols <- c(condition_col, celltype_col)",
    fixed = TRUE
  )
  expect_match(
    text,
    "Every condition-by-cell-type Pando GRN must complete successfully",
    fixed = TRUE
  )
})

test_that("Stage 3 uses active bootstrap targets rather than target projection", {
  text <- paste(
    deparse(body(.rc_build_condition_meta_modules)),
    collapse = "\n"
  )
  expect_match(text, ".rc_summarize_supported_metabolic_genes", fixed = TRUE)
  expect_match(text, "rc_map_meta_module_core_reactions", fixed = TRUE)
  expect_match(text, "bootstrap-stable", fixed = TRUE)
  expect_false(grepl("rc_project_metabolic_grn", text, fixed = TRUE))
  expect_false(grepl("top_k_neighbors", text, fixed = TRUE))
})

test_that("merged meta-modules contain biological reactions only", {
  condition_modules <- list(
    group_status = data.frame(group_id = "A|T", status = "ok"),
    tf_peak_gene_candidates = data.frame(edge_universe_id = "u1"),
    tf_peak_gene_all = data.frame(),
    tf_peak_gene_significant = data.frame(),
    supported_metabolic_genes = data.frame(),
    core_gene_reaction = data.frame(
      group_id = "A|T",
      module_id = "A|T::SUPPORTED_METABOLIC_GENES",
      reaction_id = "R1",
      is_core = TRUE
    ),
    reaction_membership = data.frame(
      group_id = "A|T",
      module_id = "A|T::SUPPORTED_METABOLIC_GENES",
      reaction_id = c("R1", "R2")
    ),
    meta_module_summary = data.frame()
  )
  out <- .rc_merge_meta_module_catalogue(condition_modules)
  expect_setequal(
    out$merged_reaction_membership$reaction_id,
    c("R1", "R2")
  )
  expect_setequal(out$merged_core_reactions$reaction_id, "R1")
  expect_false(out$is_gem)
  expect_false(out$fastcore_applied)
  expect_identical(out$source_edge_universe_ids, "u1")
  expect_identical(out$source_group_ids, "A|T")
  expect_false(any(grepl(
    "fastcore",
    out$merged_reaction_membership$inclusion_stage,
    ignore.case = TRUE
  )))
})

test_that("metacell construction is condition only and label guided", {
  text <- paste(
    deparse(body(.rc_make_condition_pooled_metacells)),
    collapse = "\n"
  )
  expect_match(
    text,
    'object@meta.data[[internal_celltype_col]] <- "all_celltypes"',
    fixed = TRUE
  )
  expect_match(text, 'pooling_scope <- "condition_only"', fixed = TRUE)
  expect_match(text, "metacell_grouping = condition_col", fixed = TRUE)
  expect_match(text, "metacell_args$gamma <- 30L", fixed = TRUE)
  expect_match(text, "label_col = celltype_col", fixed = TRUE)
  expect_false(grepl("sample_balance", text, fixed = TRUE))
  expect_false(grepl("sample_weighting", text, fixed = TRUE))
})

test_that("canonical metacells automatically use cell type as the label", {
  step_formals <- formals(rc_regcompass_step_metacells)
  run_formals <- formals(rc_run_regcompass)

  expect_false("label_col" %in% names(step_formals))
  expect_false("metacell_label_col" %in% names(run_formals))
  expect_false("sample_col" %in% names(step_formals))
  expect_false("sample_col" %in% names(run_formals))
})

test_that("dominant cell type is assigned after condition-only metacells", {
  skip_if_not_installed("SeuratObject")
  counts <- Matrix::Matrix(
    matrix(
      1,
      nrow = 1,
      ncol = 6,
      dimnames = list("g1", paste0("c", 1:6))
    ),
    sparse = TRUE
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object$cell_type <- c("T", "T", "B", "B", "B", "T")
  pooled <- list(
    membership = data.frame(
      cell_id = paste0("c", 1:6),
      metacell_id = c("m1", "m1", "m1", "m2", "m2", "m2")
    ),
    metacell_meta = data.frame(metacell_id = c("m1", "m2"))
  )
  out <- .rc_assign_metacell_dominant_celltype(
    pooled,
    object,
    "cell_type"
  )
  expect_identical(out$metacell_meta$cell_type, c("T", "B"))
  expect_equal(
    out$metacell_meta$dominant_celltype_fraction,
    c(2 / 3, 2 / 3)
  )
  expect_true(all(out$metacell_meta$mixed_celltype_metacell))
  expect_false(any(out$metacell_meta$dominant_celltype_tied))
})

test_that("condition metacells reject fragment pooling without maps", {
  text <- paste(
    deparse(body(.rc_make_condition_pooled_metacells)),
    collapse = "\n"
  )
  expect_match(
    text,
    "requires `fragment_files = FALSE`",
    fixed = TRUE
  )
})
