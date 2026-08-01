test_that("Stage 1 derives memory controls from its actual output directory", {
  outdir <- tempfile("regcompass-step1-")
  control <- .rc_pando_engine_control(
    outdir,
    list(memory_budget_mb = 768, resume = FALSE)
  )
  expect_identical(control$memory_budget_mb, 768)
  expect_identical(control$dense_max_p, 2048L)
  expect_identical(control$lambda_batch_size, 1L)
  expect_identical(control$diagnostics_level, "compact")
  expect_identical(
    control$checkpoint_dir,
    file.path(outdir, "target_checkpoints")
  )
  expect_false(control$resume)

  disabled <- .rc_pando_engine_control(
    outdir, list(checkpoint_dir = NULL)
  )
  expect_null(disabled$checkpoint_dir)

  expect_error(
    .rc_pando_engine_control(outdir, list(unknown = TRUE)),
    "Unknown"
  )
})

test_that("execution summaries preserve dense and matrix-free diagnostics", {
  diagnostics <- data.frame(
    network_name = rep("regcompass_condition_grn", 3),
    target = c("G1", "G2", "G3"),
    stage = c("complete", "complete", "target_skip"),
    path_backend = c(
      "dense_centered_gram", "sparse_matrix_free", NA_character_
    ),
    refit_backend = c(
      "dense_direct_schur", "matrix_free_schur_pcg", NA_character_
    ),
    predictors = c(100, 20000, NA),
    nonzeros = c(1000, 120000, NA),
    pcg_iterations = c(NA, 28, NA),
    pcg_residual = c(NA, 2e-9, NA),
    estimated_peak_bytes = c(1e6, 2e8, NA),
    stringsAsFactors = FALSE
  )
  summary <- .rc_pando_execution_summary(diagnostics)
  expect_identical(summary$targets_total, 3L)
  expect_identical(summary$targets_dense, 1L)
  expect_identical(summary$targets_matrix_free, 1L)
  expect_identical(summary$targets_failed, 1L)
  expect_identical(summary$largest_p, 20000L)
  expect_equal(summary$largest_nnz, 120000)
  expect_identical(summary$max_pcg_iterations, 28L)
  expect_equal(summary$max_pcg_residual, 2e-9)
  expect_equal(summary$estimated_peak_bytes, 2e8)
})

test_that("Stage 1 advertises bounded batches and fail-fast target events", {
  source <- paste(
    readLines("R/condition_grn_contract.R", warn = FALSE),
    collapse = "\n"
  )
  for (phase in c(
    "target_plan", "target_batch_start", "target_batch_complete",
    "matrix_free_target", "dense_target", "target_failure",
    "checkpoint_written"
  )) {
    expect_match(source, phase, fixed = TRUE)
  }
})
