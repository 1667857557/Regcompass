test_that("v1.8.4 public workflow is GRN first", {
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
  expect_false("inference_unit" %in% names(formals(rc_run_regcompass)))
  expect_identical(eval(formals(rc_run_regcompass)$fragment_files), FALSE)
  expect_true("meta_module_args" %in% names(formals(rc_run_regcompass)))
})

test_that("Pando grouping is condition by cell type", {
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

test_that("Stage 3 uses significant targets rather than target projection", {
  text <- paste(
    deparse(body(.rc_build_condition_meta_modules)),
    collapse = "\n"
  )
  expect_match(text, ".rc_summarize_supported_metabolic_genes", fixed = TRUE)
  expect_match(text, "rc_map_meta_module_core_reactions", fixed = TRUE)
  expect_false(grepl("rc_project_metabolic_grn", text, fixed = TRUE))
  expect_false(grepl("top_k_neighbors", text, fixed = TRUE))
})

test_that("merged meta-modules contain biological reactions only", {
  condition_modules <- list(
    sample_status = data.frame(status = "ok"),
    tf_peak_gene_all = data.frame(),
    tf_peak_gene_significant = data.frame(),
    supported_metabolic_genes = data.frame(),
    core_gene_reaction = data.frame(
      sample_id = "A|T",
      module_id = "A|T::SUPPORTED_METABOLIC_GENES",
      reaction_id = "R1",
      is_core = TRUE
    ),
    reaction_membership = data.frame(
      sample_id = "A|T",
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
  expect_false(any(grepl(
    "fastcore",
    out$merged_reaction_membership$inclusion_stage,
    ignore.case = TRUE
  )))
})

test_that("metacell construction is condition-only without sample balancing", {
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
  expect_match(text, 'sample_weighting <- "none"', fixed = TRUE)
  expect_match(text, "metacell_grouping = condition_col", fixed = TRUE)
  expect_match(text, "gamma <- 30L", fixed = TRUE)
  expect_match(text, "Sample balancing is not part", fixed = TRUE)
  expect_match(text, "label_col = celltype_col", fixed = TRUE)
  expect_match(text, '"label_col"', fixed = TRUE)
  expect_false(grepl("label_col = label_col", text, fixed = TRUE))
})

test_that("canonical metacells automatically use cell type as the label", {
  step_formals <- formals(rc_regcompass_step_metacells)
  run_formals <- formals(rc_run_regcompass)

  expect_false("label_col" %in% names(step_formals))
  expect_false("metacell_label_col" %in% names(run_formals))
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
