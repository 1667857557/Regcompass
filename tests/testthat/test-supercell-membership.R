test_that("metacell metadata may mix biological samples but not biology strata", {
  membership <- data.frame(
    cell_id = c("c1", "c2"),
    metacell_id = c("MC1", "MC1"),
    sample_id = c("S1", "S2"),
    condition = c("Control", "Control"),
    cell_type = c("T", "T"),
    stringsAsFactors = FALSE
  )
  out <- rc_build_metacell_metadata(membership)
  expect_equal(out$n_cells, 2L)
  expect_error(
    rc_build_metacell_metadata(transform(
      membership, cell_type = c("T", "B")
    )),
    "cell_type"
  )
})
