test_that("public API contains only the supported workflow", {
  expect_setequal(
    getNamespaceExports("RegCompassR"),
    c(
      "rc_prepare_human2_gem",
      "rc_make_medium_scenarios",
      "rc_run_regcompass",
      "rc_run_regcompass_one_shot"
    )
  )
})


test_that("deprecated and versioned entry points are absent", {
  retired <- c(
    "rc_prepare_human2_gem_v12",
    "rc_run_regcompass_multiome_metacell",
    "rc_run_layer1_from_metacells",
    "rc_recompute_metacell_peak_gene_links",
    "rc_recompute_metacell_peak_gene_links_by_stratum",
    "rc_build_meta_module_gem_cache"
  )
  expect_false(any(vapply(retired, exists, logical(1), inherits = TRUE)))

  source_dir <- if (dir.exists("R")) {
    "R"
  } else {
    normalizePath(file.path("..", "..", "R"), mustWork = TRUE)
  }
  legacy_late_files <- basename(list.files(source_dir, pattern = "^zzz"))
  expect_length(legacy_late_files, 0L)

  workflow_stages <- c(
    "workflow_stage_01_architecture.R",
    "workflow_stage_02_compatibility.R",
    "workflow_stage_03_signed_projection.R",
    "workflow_stage_04_result_contracts.R",
    "workflow_stage_04b_human_medium_presets.R",
    "workflow_stage_05_api_contracts.R",
    "workflow_stage_06_audit_contracts.R"
  )
  expect_true(all(file.exists(file.path(source_dir, workflow_stages))))

  description <- read.dcf(file.path(source_dir, "..", "DESCRIPTION"))
  collate <- strsplit(description[1L, "Collate"], "[[:space:]]+")[[1L]]
  collate <- gsub("^['\"]|['\"]$", "", collate[nzchar(collate)])
  expect_identical(utils::tail(collate, length(workflow_stages)), workflow_stages)
})
