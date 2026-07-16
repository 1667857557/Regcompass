test_that("strict-stratum metadata are validated", {
  expect_identical(
    .rc_strict_stratum_cols("sample_id", "condition", "cell_type"),
    c("condition", "sample_id", "cell_type")
  )
  expect_error(
    .rc_strict_stratum_cols("sample_id", "sample_id", "cell_type"),
    "distinct"
  )
  meta <- data.frame(
    condition = c("ctrl", "ctrl"),
    sample_id = c("s1", ""),
    cell_type = c("T", "T")
  )
  expect_error(
    .rc_add_stratum_id(meta, c("condition", "sample_id", "cell_type")),
    "missing"
  )
})

test_that("stratum IDs use the shared pipe-separated convention", {
  meta <- data.frame(
    condition = c("ctrl", "stim"),
    sample_id = c("s1", "s2"),
    cell_type = c("T", "B")
  )
  expect_identical(
    rc_make_stratum_id(meta, c("condition", "sample_id", "cell_type")),
    c("ctrl|s1|T", "stim|s2|B")
  )
})

test_that("Pando uses strict-stratum grouping", {
  body_text <- paste(deparse(body(rc_run_pando_meta_modules)), collapse = "\n")
  expect_match(body_text, "rc_make_stratum_id", fixed = TRUE)
  expect_match(body_text, "rc_parallel_lapply(group_ids", fixed = TRUE)
})
