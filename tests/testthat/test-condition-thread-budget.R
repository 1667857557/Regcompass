test_that("condition GRN uses one parallel level at a time", {
  body_text <- paste(
    deparse(body(.rc_fit_condition_grns_by_cell_type)), collapse = "\n"
  )
  expect_match(
    body_text,
    "inner_parallel <- condition_parallel && length(fit_tasks) == 1L",
    fixed = TRUE
  )

  fit_text <- paste(
    deparse(body(.rc_condition_multitask_fit_task)), collapse = "\n"
  )
  expect_false(grepl("REGCOMPASS_TASK_WORKERS", fit_text, fixed = TRUE))
  expect_false(grepl("nested_pool", fit_text, fixed = TRUE))
  expect_false(exists(
    ".rc_condition_nested_target_bpparam", inherits = TRUE
  ))
  expect_false(exists(
    ".rc_allocate_condition_target_workers", inherits = TRUE
  ))
})

test_that("generic parallel dispatcher does not capture full task lists", {
  body_text <- paste(deparse(body(rc_parallel_lapply)), collapse = "\n")
  expect_false(grepl("REGCOMPASS_DISPATCH", body_text, fixed = TRUE))
  expect_false(grepl("snow_hierarchical", body_text, fixed = TRUE))
  expect_false(grepl("task_keys", body_text, fixed = TRUE))
  expect_false(grepl("budgets", body_text, fixed = TRUE))
})

test_that("Stage 1 Pando dispatch trims unused Seurat payload", {
  helper_text <- paste(
    deparse(body(.rc_stage1_pando_working_object)), collapse = "\n"
  )
  expect_match(helper_text, "Seurat::DietSeurat", fixed = TRUE)
  expect_match(helper_text, "counts = TRUE", fixed = TRUE)
  expect_match(helper_text, "data = TRUE", fixed = TRUE)
  expect_match(helper_text, "scale.data = FALSE", fixed = TRUE)
  expect_match(helper_text, "dimreducs = NULL", fixed = TRUE)
  expect_match(helper_text, "graphs = NULL", fixed = TRUE)
  expect_match(helper_text, "misc = FALSE", fixed = TRUE)

  route_text <- paste(
    deparse(body(.rc_run_condition_pando_batch)), collapse = "\n"
  )
  expect_match(
    route_text, ".rc_stage1_pando_working_object", fixed = TRUE
  )
  expect_false(grepl("condition_object <- subset", route_text, fixed = TRUE))
})

test_that("condition workers release superseded large objects", {
  prepare_text <- paste(
    deparse(body(.rc_condition_prepare_celltype_task)), collapse = "\n"
  )
  expect_match(prepare_text, "task$object <- NULL", fixed = TRUE)
  expect_match(prepare_text, "one <- NULL", fixed = TRUE)
  expect_match(prepare_text, "motif <- NULL", fixed = TRUE)

  fit_text <- paste(
    deparse(body(.rc_condition_multitask_fit_task)), collapse = "\n"
  )
  expect_match(fit_text, "task$grn <- NULL", fixed = TRUE)
  expect_match(fit_text, "args <- NULL", fixed = TRUE)
})
