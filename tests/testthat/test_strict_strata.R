test_that("public stratum metadata are validated", {
  expect_error(
    rc_make_stratum_id(data.frame(), character()),
    "at least one metadata column"
  )
  meta <- data.frame(
    condition = c("ctrl", "ctrl"),
    sample_id = c("s1", ""),
    cell_type = c("T", "T")
  )
  expect_error(
    rc_make_stratum_id(meta, c("condition", "sample_id", "cell_type")),
    "missing or empty"
  )
  expect_error(
    rc_make_stratum_id(meta, c("condition", "missing")),
    "Missing stratum columns"
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
