test_that("condition scheduler plans global plus condition x cell-type tasks", {
  metadata <- data.frame(
    condition = rep(c("A", "B"), each = 6L, times = 2L),
    cell_type = rep(c("T", "B"), each = 12L),
    row.names = paste0("cell", seq_len(24L)),
    stringsAsFactors = FALSE
  )
  plan <- .rc_condition_parallel_plan(
    metadata = metadata,
    condition_types = c("T", "B"),
    condition_col = "condition",
    celltype_col = "cell_type",
    min_cells = 3L
  )
  expect_identical(names(plan), c("T", "B"))
  expect_identical(plan$T$conditions, c("A", "B"))
  expect_identical(plan$B$conditions, c("A", "B"))
  expect_length(plan$T$global_cells, 12L)
  expect_length(plan$T$cells_by_condition$A, 6L)
  expect_length(plan$T$cells_by_condition$B, 6L)
})

test_that("condition scheduler preserves the min-cells error contract", {
  metadata <- data.frame(
    condition = c(rep("A", 5L), rep("B", 2L)),
    cell_type = "T",
    row.names = paste0("cell", seq_len(7L)),
    stringsAsFactors = FALSE
  )
  expect_error(
    .rc_condition_parallel_plan(
      metadata = metadata,
      condition_types = "T",
      condition_col = "condition",
      celltype_col = "cell_type",
      min_cells = 3L
    ),
    "below min_cells"
  )
})

test_that("condition penalty eligibility requires an ok target fit plus estimable BH significance", {
  coefficient <- data.frame(
    estimate = c(1e-6, 2, 3, 4),
    padj = c(0.01, 0.01, 0.051, 0.01),
    estimable = c(TRUE, TRUE, TRUE, TRUE),
    fit_status = c("ok", "rank_deficient", "ok", "insufficient_df"),
    corr = c(0, 1, 1, 1),
    stringsAsFactors = FALSE
  )
  expect_identical(
    .rc_condition_penalty_gate(coefficient),
    c(TRUE, FALSE, FALSE, FALSE)
  )
})

test_that("condition fit status maps exactly by target and condition", {
  coefficient <- data.frame(
    target = c("G1", "G1", "G2", "G2"),
    condition = c("A", "B", "A", "B"),
    stringsAsFactors = FALSE
  )
  fit <- list(
    fit = data.frame(
      target = c("G2", "G1", "G2", "G1"),
      condition = c("B", "A", "A", "B"),
      fit_status = c("rank_deficient", "ok", "ok", "insufficient_df"),
      stringsAsFactors = FALSE
    )
  )
  expect_identical(
    .rc_condition_fit_status_for_coefficients(fit, coefficient),
    c("ok", "insufficient_df", "ok", "rank_deficient")
  )
})

test_that("standard Pando post-fit filter is BH-only after candidate screening", {
  table <- data.frame(
    estimate = c(1e-8, 1, 1),
    padj = c(0.01, 0.049, 0.051),
    corr = c(0, 1, 1),
    estimable = c(TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  observed <- .rc_filter_standard_pando_edges(table)
  expect_equal(nrow(observed), 2L)
  expect_equal(observed$padj, c(0.01, 0.049))
})

test_that("RegCompass condition scheduling uses exported Pando primitives", {
  skip_if_not_installed("Pando")
  exports <- getNamespaceExports("Pando")
  expect_true(all(c(
    "discover_grn_edges", "union_grn_edges", "fit_grn_from_edges"
  ) %in% exports))
})

test_that("condition route calls the canonical RegCompass condition function", {
  body_text <- paste(deparse(body(.rc_run_condition_pando_batch)), collapse = "\n")
  expect_match(body_text, "rc_fit_condition_grns_by_cell_type")
  expect_false(grepl("parallel_scope", body_text, fixed = TRUE))
  expect_false(grepl("rc_fit_condition_grns_regcompass_parallel", body_text,
                     fixed = TRUE))
})

test_that("Pando workflow has no versioned or wrapper implementation functions", {
  forbidden <- c(
    ".rc_route_pando_infer_args_core",
    ".rc_extract_condition_grn_contract_core",
    ".rc_merge_condition_job_results_core",
    ".rc_merge_pando_results_core",
    ".rc_merge_pando_results_validated",
    ".rc_merge_pando_results_with_parallel_objects",
    ".rc_fit_condition_grns_regcompass_parallel"
  )
  expect_false(any(vapply(forbidden, exists, logical(1), inherits = TRUE)))
})
