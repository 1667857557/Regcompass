test_that("v2.0.0 public workflow is GRN first", {
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

test_that("Pando shares one cell-type fit across conditions", {
  implementation <- paste(
    deparse(body(.rc_fit_condition_grns_by_cell_type)), collapse = "\n"
  )
  expect_match(
    implementation,
    "group_cols <- c(condition_col, celltype_col)",
    fixed = TRUE
  )
  expect_match(implementation, "Pando::infer_condition_grn", fixed = TRUE)
  expect_match(implementation, "condition_grn_fits", fixed = TRUE)
  expect_match(implementation, "cell_type = cell_type", fixed = TRUE)
})

test_that("Stage 3 uses active targets rather than target projection", {
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
    condition_fit_status = data.frame(status = "ok"),
    tf_peak_gene_condition_all = data.frame(),
    tf_peak_gene_condition = data.frame(),
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

test_that("metacells use only condition and broad cell type as hard strata", {
  text <- paste(
    deparse(body(.rc_make_condition_celltype_metacells)),
    collapse = "\n"
  )
  expect_match(text, 'pooling_scope <- "condition_by_cell_type"', fixed = TRUE)
  expect_match(
    text, "metacell_grouping = c(condition_col, celltype_col)",
    fixed = TRUE
  )
  expect_match(text, "gamma <- 30L", fixed = TRUE)
  expect_match(
    text, "pooled$membership[[supercell_stratum_col]] <- NULL", fixed = TRUE
  )
  expect_match(text, "label_col = NULL", fixed = TRUE)
  expect_match(text, '"label_col"', fixed = TRUE)
  expect_false(grepl("label_col = label_col", text, fixed = TRUE))
})

test_that("canonical metacells automatically use cell type as the label", {
  step_formals <- formals(rc_regcompass_step_metacells)
  run_formals <- formals(rc_run_regcompass)

  expect_false("label_col" %in% names(step_formals))
  expect_false("metacell_label_col" %in% names(run_formals))
})

test_that("canonical construction does not use posthoc cell-type assignment", {
  text <- paste(
    deparse(body(.rc_make_condition_celltype_metacells)),
    collapse = "\n"
  )
  expect_false(grepl(
    ".rc_assign_metacell_dominant_celltype", text, fixed = TRUE
  ))
  expect_match(
    text, "mixing condition or broad cell type", fixed = TRUE
  )
})

test_that("condition metacells reject fragment pooling without maps", {
  text <- paste(
    deparse(body(.rc_make_condition_celltype_metacells)),
    collapse = "\n"
  )
  expect_match(
    text,
    "requires `fragment_files = FALSE`",
    fixed = TRUE
  )
})
