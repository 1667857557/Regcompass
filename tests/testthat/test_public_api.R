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
  late_files <- basename(list.files(source_dir, pattern = "^zzz"))
  expect_setequal(
    late_files,
    c(
      "zzz_architecture_correctness.R",
      "zzzz_architecture_hotfixes.R",
      "zzzzz_signed_projection.R",
      "zzzzzz_signal_scale.R",
      "zzzzzzzz_required_result_fixes.R",
      "zzzzzzzzz_required_medium_hotfix.R",
      "zzzzzzzzzz_required_contract_hotfix.R",
      "zzzzzzzzzzz_force_medium_application.R"
    )
  )
})
