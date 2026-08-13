test_that("Step 2 model batch scheduler balances unequal target counts", {
  rows_a <- paste0("A", seq_len(3223L))
  rows_b <- paste0("B", seq_len(4484L))
  model_keys <- c(
    stats::setNames(rep("model_A", length(rows_a)), rows_a),
    stats::setNames(rep("model_B", length(rows_b)), rows_b)
  )
  tasks <- RegCompassR:::.rc_step2_model_batches(
    model_keys, workers = 80L
  )
  expect_length(tasks, 80L)
  observed <- unlist(lapply(tasks, `[[`, "row_ids"), use.names = FALSE)
  expect_identical(length(observed), length(unique(observed)))
  expect_setequal(observed, names(model_keys))

  sizes <- vapply(tasks, function(task) length(task$row_ids), integer(1))
  expect_lte(max(sizes) - min(sizes), 4L)
  task_models <- vapply(tasks, `[[`, character(1), "model_key")
  expect_gt(sum(task_models == "model_B"), sum(task_models == "model_A"))
})

test_that("Step 2 model batch scheduler keeps one batch per model when workers are scarce", {
  model_keys <- stats::setNames(
    c("m1", "m2", "m3", "m4"), paste0("row", seq_len(4L))
  )
  tasks <- RegCompassR:::.rc_step2_model_batches(
    model_keys, workers = 2L
  )
  expect_length(tasks, 4L)
  expect_setequal(
    unlist(lapply(tasks, `[[`, "row_ids"), use.names = FALSE),
    names(model_keys)
  )
})
