test_that("stage progress reports percent, elapsed time and structured context", {
  outdir <- tempfile("regcompass-progress-")
  monitor <- .rc_step_monitor_start(
    "grn", outdir, progress = TRUE, total_parts = 4L
  )
  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)

  messages <- testthat::capture_messages(
    .rc_step_monitor_event(
      monitor,
      phase = "nested_cv",
      detail = "running fused target engine",
      current = 2L,
      context = list(
        cell_types = 3L,
        conditions = 2L,
        targets = 120L,
        outer_folds = 5L,
        inner_folds = 5L,
        nlambda = 50L
      )
    )
  )
  text <- paste(messages, collapse = "\n")
  expect_match(text, "phase=nested_cv", fixed = TRUE)
  expect_match(text, "50.0%", fixed = TRUE)
  expect_match(text, "elapsed=", fixed = TRUE)
  expect_match(text, "targets=120", fixed = TRUE)

  log <- utils::read.delim(
    file.path(outdir, "step_progress.tsv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  expect_true(all(c(
    "sequence", "timestamp", "stage", "phase", "status", "current",
    "total", "percent", "elapsed_seconds", "elapsed_hms", "detail",
    "context"
  ) %in% colnames(log)))
  expect_true(all(c("stage_start", "nested_cv") %in% log$phase))
  expect_match(log$context[log$phase == "nested_cv"], "outer_folds=5")
})

test_that("progress FALSE suppresses messages but retains the audit log", {
  outdir <- tempfile("regcompass-progress-silent-")
  expect_silent({
    monitor <- .rc_step_monitor_start(
      "grn", outdir, progress = FALSE, total_parts = 3L
    )
    .rc_step_monitor_event(
      monitor, "motif_mapping", "mapped motifs", current = 2L,
      context = list(targets = 10L)
    )
    .rc_step_monitor_fail(monitor)
  })
  log <- utils::read.delim(
    file.path(outdir, "step_progress.tsv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  expect_true(all(c("stage_start", "motif_mapping", "stage_error") %in%
    log$phase))
})

test_that("Stage 1 source forwards monitor and verbose status to Pando", {
  source_text <- paste(
    readLines(test_path("..", "..", "R", "stepwise_workflow.R"), warn = FALSE),
    collapse = "\n"
  )
  bridge_text <- paste(
    readLines(test_path("..", "..", "R", "condition_grn_contract.R"),
      warn = FALSE),
    collapse = "\n"
  )
  expect_match(source_text, "call_args\\$progress_monitor <- monitor")
  expect_match(source_text, "infer_args\\$verbose <- infer_args\\$verbose")
  expect_match(bridge_text, 'progress_monitor, "nested_cv"', fixed = TRUE)
  expect_match(bridge_text, "exact_sufficient_statistics")
  expect_match(bridge_text, "outer_selected_model_only")
})
