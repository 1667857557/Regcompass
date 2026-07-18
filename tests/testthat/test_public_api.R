test_that("public API contains the canonical species-aware workflow", {
  expect_setequal(
    getNamespaceExports("RegCompassR"),
    c(
      "rc_prepare_gem",
      "rc_prepare_human2_gem",
      "rc_prepare_mouse_gem",
      "rc_make_medium_scenarios",
      "rc_run_regcompass",
      "rc_run_regcompass_one_shot"
    )
  )
})


test_that("staged override architecture has been removed", {
  description <- utils::packageDescription("RegCompassR")
  expect_false(grepl(
    "workflow_stage_", description$Collate %||% "", fixed = TRUE
  ))

  candidates <- c("R", file.path("..", "R"), file.path("..", "..", "R"))
  candidates <- candidates[dir.exists(candidates)]
  if (!length(candidates)) {
    skip("Source R directory is unavailable in the installed-package test context.")
  }
  source_dir <- normalizePath(candidates[[1L]], mustWork = TRUE)
  expect_length(
    list.files(source_dir, pattern = "^workflow_stage_", full.names = TRUE),
    0L
  )

  canonical <- c(
    "rc_prepare_gem", "rc_prepare_human2_gem", "rc_prepare_mouse_gem",
    "rc_make_medium_scenarios", "rc_apply_medium_constraints",
    "rc_q95_calibrate", "rc_reaction_capacity", "rc_run_microcompass",
    "rc_run_regcompass", "rc_run_regcompass_one_shot"
  )
  source_text <- vapply(
    list.files(source_dir, pattern = "[.]R$", full.names = TRUE),
    function(path) paste(readLines(path, warn = FALSE), collapse = "\n"),
    character(1)
  )
  for (name in canonical) {
    pattern <- paste0("(?m)^", gsub("[.]", "\\\\.", name),
                      "[[:space:]]*<-[[:space:]]*function[[:space:]]*\\(")
    expect_equal(
      sum(vapply(source_text, function(x) {
        matches <- gregexpr(pattern, x, perl = TRUE)[[1L]]
        any(matches > 0L)
      }, logical(1))),
      1L,
      info = paste("canonical function", name, "must have one source definition")
    )
  }
})


test_that("retired entry points remain absent", {
  retired <- c(
    "rc_prepare_human2_gem_v12",
    "rc_run_regcompass_multiome_metacell",
    "rc_run_layer1_from_metacells",
    "rc_recompute_metacell_peak_gene_links",
    "rc_recompute_metacell_peak_gene_links_by_stratum",
    "rc_build_meta_module_gem_cache"
  )
  expect_false(any(vapply(retired, exists, logical(1), inherits = TRUE)))
})
