test_that("condition scheduler plans pooled/global plus condition cell sets", {
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

test_that("candidate provenance does not veto a BH-active Pando edge", {
  estimate <- c(1e-6, 2, 3, 4)
  padj <- c(0.01, 0.01, 0.051, 0.01)
  estimable <- rep(TRUE, 4L)
  statistically_supported <- estimable & is.finite(padj) & padj < 0.05
  global_support <- c(TRUE, TRUE, TRUE, FALSE)
  local_support <- c(FALSE, FALSE, TRUE, FALSE)
  active <- statistically_supported
  coefficient <- data.frame(
    estimate = estimate,
    padj = padj,
    estimable = estimable,
    statistically_supported = statistically_supported,
    global_support = global_support,
    local_support = local_support,
    active = active,
    significant = active,
    penalty_effect = ifelse(active, estimate, 0),
    fit_status = c("ok", "rank_deficient", "ok", "ok"),
    rsq = c(0.8, 0.9, 0.9, 0.8),
    stringsAsFactors = FALSE
  )
  expect_identical(
    .rc_condition_penalty_gate(coefficient),
    c(TRUE, FALSE, FALSE, TRUE)
  )
})

test_that("condition fit diagnostics map exactly by target and condition", {
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
      rsq = c(0.2, 0.9, 0.7, NA_real_),
      stringsAsFactors = FALSE
    )
  )
  diagnostics <- .rc_condition_fit_diagnostics_for_coefficients(fit, coefficient)
  expect_identical(
    diagnostics$fit_status,
    c("ok", "insufficient_df", "ok", "rank_deficient")
  )
  expect_equal(diagnostics$target_rsq, c(0.9, NA, 0.7, 0.2))
})

test_that("standard Pando post-fit filter uses BH plus target-model R2", {
  table <- data.frame(
    estimate = c(1e-8, 1, 1, 2),
    padj = c(0.01, 0.049, 0.051, 0.01),
    rsq = c(0.8, 0.2, 0.9, 0.01),
    corr = c(0, 1, 1, 1),
    estimable = c(TRUE, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  observed <- .rc_filter_standard_pando_edges(
    table, target_rsq_threshold = 0.05
  )
  expect_equal(nrow(observed), 2L)
  expect_equal(observed$padj, c(0.01, 0.049))
})

test_that("RegCompass condition scheduling uses Pando condition GRN API", {
  skip_if_not_installed("Pando")
  exports <- getNamespaceExports("Pando")
  expect_true("infer_condition_grn" %in% exports)
  expect_true("condition_grn_fit" %in% exports)
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
    ".rc_fit_condition_grns_regcompass_parallel",
    ".rc_condition_multitask_fit_task"
  )
  expect_false(any(vapply(forbidden, exists, logical(1), inherits = TRUE)))
})
