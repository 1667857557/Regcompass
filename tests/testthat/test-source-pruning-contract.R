test_that("retired source subsystems stay out of the package", {
  description <- utils::packageDescription("RegCompassR")
  collate <- description$Collate %||% ""
  retired <- c(
    "metacell_import.R",
    "metacell_peak_calling.R",
    "metacell_object_merge.R",
    "metacell_fragments.R",
    "condition_direction_report.R"
  )
  expect_false(any(vapply(
    retired, grepl, logical(1), x = collate, fixed = TRUE
  )))
  expect_true(grepl("seurat_fragments.R", collate, fixed = TRUE))
  expect_false("rc_report_condition_directions" %in%
                 getNamespaceExports("RegCompassR"))
})

test_that("retired workflow controls stay removed", {
  expect_false("fragment_files" %in% names(formals(rc_run_regcompass)))
  expect_false("fragment_files" %in%
                 names(formals(rc_run_regcompass_one_shot)))
  expect_false("fragment_files" %in%
                 names(formals(rc_regcompass_step_grn)))
  expect_false("fragment_files" %in%
                 names(formals(rc_regcompass_step_metacells)))
  expect_false("fragment_files" %in%
                 names(formals(.rc_make_condition_celltype_metacells)))
  expect_false("overwrite" %in% names(.rc_condition_metacell_defaults()))
})

test_that("current core validation and mapping helpers remain", {
  ns <- asNamespace("RegCompassR")
  expect_true(exists(".rc_validate_seurat_stack_versions", envir = ns,
                     inherits = FALSE))
  expect_true(exists("rc_validate_gem", envir = ns, inherits = FALSE))
  expect_true(exists("rc_align_bound", envir = ns, inherits = FALSE))
  expect_true(exists(".rc_medium_exchange_metabolites", envir = ns,
                     inherits = FALSE))
  expect_true(exists(".rc_read_gem", envir = ns, inherits = FALSE))
  expect_false(exists("rc_make_gem", envir = ns, inherits = FALSE))
  expect_false(exists(".rc_medium_resolve_exchange_metabolites", envir = ns,
                      inherits = FALSE))
})
