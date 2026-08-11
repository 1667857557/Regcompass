test_that("Pando routing defaults peak_cor to 0.05", {
  routed <- .rc_route_pando_infer_args(
    list(),
    condition_types = character(),
    standard_types = "T_cell"
  )

  expect_equal(routed$standard$tf_cor, 0.1)
  expect_equal(routed$standard$peak_cor, 0.05)
  expect_equal(routed$standard$adjust_method, "BH")
})

test_that("single-condition standard Pando drops condition-only controls", {
  routed <- .rc_route_pando_infer_args(
    list(
      tf_cor = 0.1,
      peak_cor = 0,
      adjust_method = "BH",
      padj_threshold = 0.05,
      rank_action = "mark",
      min_residual_df = 1L
    ),
    condition_types = character(),
    standard_types = "T_cell"
  )

  expect_equal(routed$standard$tf_cor, 0.1)
  expect_equal(routed$standard$peak_cor, 0)
  expect_equal(routed$standard$adjust_method, "BH")
  expect_false(any(c(
    "padj_threshold", "rank_action", "min_residual_df"
  ) %in% names(routed$standard)))
  expect_setequal(
    routed$diagnostics$argument,
    c("padj_threshold", "rank_action", "min_residual_df")
  )
})

test_that("condition GRN routes ridge controls and drops standard-model controls", {
  ridge <- list(
    lambda_grid = c(0.01, 0.1, 1),
    lambda_rule = "1se",
    fusion_ratio = 1,
    cv_folds = 4L,
    seed = 7L,
    scale_floor = 1e-8
  )
  routed <- .rc_route_pando_infer_args(
    list(
      tf_cor = 0.2,
      peak_cor = 0.05,
      adjust_method = "BH",
      padj_threshold = 0.2,
      rank_action = "mark",
      min_residual_df = 2L,
      condition_ridge_control = ridge,
      method = "glmnet",
      alpha = 0.5,
      scale = TRUE
    ),
    condition_types = "Monocyte",
    standard_types = character()
  )

  expect_equal(routed$condition$tf_cor, 0.2)
  expect_equal(routed$condition$peak_cor, 0.05)
  expect_equal(routed$condition$padj_threshold, 0.2)
  expect_equal(routed$condition$min_residual_df, 2L)
  expect_identical(routed$condition$condition_ridge_control, ridge)
  expect_false(any(c("method", "alpha", "scale") %in%
                     names(routed$condition)))
  expect_setequal(
    routed$diagnostics$argument,
    c("method", "alpha", "scale")
  )
})

test_that("mixed routing preserves controls for their own route", {
  routed <- .rc_route_pando_infer_args(
    list(
      tf_cor = 0.1,
      peak_cor = 0,
      adjust_method = "BH",
      padj_threshold = 0.05,
      rank_action = "mark",
      min_residual_df = 1L,
      condition_ridge_control = list(fusion_ratio = 2),
      method = "glm",
      scale = FALSE
    ),
    condition_types = "T_cell",
    standard_types = "HSPC"
  )

  expect_true(all(c(
    "padj_threshold", "rank_action", "min_residual_df",
    "condition_ridge_control"
  ) %in% names(routed$condition)))
  expect_true(all(c("method", "scale") %in% names(routed$standard)))
  expect_false("method" %in% names(routed$condition))
  expect_false("rank_action" %in% names(routed$standard))
  expect_false("condition_ridge_control" %in% names(routed$standard))
})

test_that("unknown Pando arguments fail before model fitting", {
  expect_error(
    .rc_route_pando_infer_args(
      list(tf_cor = 0.1, misspelled_argument = TRUE),
      standard_types = "T_cell"
    ),
    "Unsupported `pando_infer_args`"
  )
})

test_that("condition BH threshold defaults to 0.05 and accepts any value in (0,1)", {
  default <- .rc_route_pando_infer_args(
    list(adjust_method = "BH"),
    condition_types = "T_cell"
  )
  expect_equal(default$condition$padj_threshold, 0.05)

  routed <- .rc_route_pando_infer_args(
    list(adjust_method = "BH", padj_threshold = 0.5),
    condition_types = "T_cell"
  )
  expect_equal(routed$condition$padj_threshold, 0.5)

  expect_error(
    .rc_route_pando_infer_args(
      list(adjust_method = "BH", padj_threshold = 1),
      condition_types = "T_cell"
    ),
    "padj_threshold in \\(0, 1\\)"
  )
})

test_that("parallel condition jobs preserve separate Pando objects and contrasts", {
  make_job <- function(cell_type, cell_id) {
    fit <- list(cell_type = cell_type)
    list(
      pando_grn_data = structure(list(id = cell_type), class = "GRNData"),
      condition_grn_fits = stats::setNames(list(fit), cell_type),
      condition_fit_status = data.frame(cell_type = cell_type),
      pando_network_index = data.frame(cell_type = cell_type),
      pando_fit_diagnostics = data.frame(cell_type = cell_type),
      tf_peak_gene_universal = data.frame(cell_type = cell_type),
      tf_peak_gene_condition_all = data.frame(cell_type = cell_type),
      tf_peak_gene_condition = data.frame(cell_type = cell_type),
      tf_peak_gene_condition_effect_all = data.frame(cell_type = cell_type),
      tf_peak_gene_condition_effect = data.frame(cell_type = cell_type),
      tf_peak_gene_condition_contrasts = data.frame(
        cell_type = cell_type, contrast = "A-B"
      ),
      paired_cell_metadata = data.frame(
        cell_id = cell_id,
        condition = "Control",
        cell_type = cell_type,
        stringsAsFactors = FALSE
      ),
      paired_cell_ids = cell_id,
      target_metabolic_genes = paste0("GENE_", cell_type),
      pando_execution_summary = list(
        fit_engine = "two_stage_exact_edge_union_multitask_ridge",
        targets_total = 1L,
        targets_failed = 0L
      )
    )
  }

  jobs <- list(
    make_job("T_cell", "T1"),
    make_job("Monocyte", "M1")
  )
  merged <- .rc_merge_condition_job_results(jobs)

  expect_setequal(
    names(merged$pando_grn_data_by_cell_type),
    c("T_cell", "Monocyte")
  )
  expect_null(merged$pando_grn_data)
  expect_setequal(merged$paired_cell_ids, c("T1", "M1"))
  expect_setequal(
    names(merged$condition_grn_fits),
    c("T_cell", "Monocyte")
  )
  expect_equal(nrow(merged$tf_peak_gene_condition_contrasts), 2L)
})
