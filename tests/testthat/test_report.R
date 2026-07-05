
test_that("report includes GPR-aware confidence and all-missing Q95 diagnostics", {
  f <- tempfile(fileext = ".md")
  q95 <- data.frame(reaction_id = "R1", all_missing_reaction_flag = TRUE, q95_power_class = "very_low")
  conf <- data.frame(
    reaction_id = "R1",
    pool_id = "p1",
    reaction_confidence = NA_real_,
    reaction_confidence_method = "gpr_aware",
    confidence_source = "gpr_aware_rna_detection",
    no_complete_gpr_group_flag = TRUE,
    reaction_unsupported_by_complete_gpr_flag = TRUE,
    any_incomplete_gpr_group_flag = TRUE,
    complete_and_group_fraction = 0,
    best_and_group_observed_fraction = 0,
    stringsAsFactors = FALSE
  )
  rc_write_report_md(f, q95_diagnostics = q95, confidence = conf)
  txt <- paste(readLines(f), collapse = "\n")
  expect_true(grepl("reaction_confidence_method", txt))
  expect_true(grepl("confidence_source", txt))
  expect_true(grepl("no_complete_gpr_group_flag", txt))
  expect_true(grepl("Complete AND-group fraction", txt))
  expect_true(grepl("all_missing_reaction_flag", txt))
})
