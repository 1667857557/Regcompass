test_that("public API contains only the canonical species-aware workflow", {
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

test_that("canonical source architecture has one definition per function", {
  description <- utils::packageDescription("RegCompassR")
  collate <- description$Collate %||% ""
  expect_false(grepl("workflow_stage_", collate, fixed = TRUE))
  expect_false(grepl("zzz", collate, fixed = TRUE))
  expect_false(grepl("calibration_q95.R", collate, fixed = TRUE))
  expect_false(grepl("global_workflow.R", collate, fixed = TRUE))
  expect_false(grepl("pando_confidence.R", collate, fixed = TRUE))
  expect_match(collate, "workflow_utils.R", fixed = TRUE)
  expect_match(collate, "pando_evidence_utils.R", fixed = TRUE)

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
  expect_length(
    list.files(source_dir, pattern = "^z+.*[.]R$", full.names = TRUE),
    0L
  )
  expect_length(
    list.files(
      source_dir,
      pattern = "^(calibration_q95|global_workflow|pando_confidence)[.]R$",
      full.names = TRUE
    ),
    0L
  )

  canonical <- c(
    "rc_prepare_gem", "rc_prepare_human2_gem", "rc_prepare_mouse_gem",
    "rc_make_medium_scenarios", "rc_apply_medium_constraints",
    "rc_reaction_capacity", "rc_compute_multiome_penalty",
    "rc_run_microcompass", "rc_run_regcompass",
    "rc_run_regcompass_one_shot"
  )
  source_text <- vapply(
    list.files(source_dir, pattern = "[.]R$", full.names = TRUE),
    function(path) paste(readLines(path, warn = FALSE), collapse = "\n"),
    character(1)
  )
  for (name in canonical) {
    pattern <- paste0(
      "(?m)^", gsub("[.]", "\\\\.", name),
      "[[:space:]]*<-[[:space:]]*function[[:space:]]*\\("
    )
    expect_equal(
      sum(vapply(source_text, function(x) {
        matches <- gregexpr(pattern, x, perl = TRUE)[[1L]]
        any(matches > 0L)
      }, logical(1))),
      1L,
      info = paste("canonical function", name, "must have one source definition")
    )
  }

  source_files <- list.files(source_dir, pattern = "[.]R$", full.names = TRUE)
  definitions <- unlist(lapply(source_files, function(path) {
    lines <- readLines(path, warn = FALSE)
    lines <- lines[grepl(
      "^[A-Za-z.][A-Za-z0-9._]*[[:space:]]*<-[[:space:]]*function",
      lines
    )]
    sub(
      "^([A-Za-z.][A-Za-z0-9._]*)[[:space:]]*<-.*$",
      "\\1",
      lines
    )
  }), use.names = FALSE)
  expect_false(
    anyDuplicated(definitions) > 0L,
    info = "every top-level function must have exactly one source definition"
  )
})

test_that("retired entry points and evidence APIs remain absent", {
  retired <- c(
    "rc_prepare_human2_gem_v12",
    "rc_download_humangem_gpr_table",
    "rc_prepare_humangem_gpr_table",
    "rc_run_regcompass_multiome_metacell",
    "rc_run_layer1_from_metacells",
    "rc_recompute_metacell_peak_gene_links",
    "rc_recompute_metacell_peak_gene_links_by_stratum",
    "rc_build_meta_module_gem_cache",
    "rc_replace_humangem_gene_ids",
    "rc_convert_humangem_yaml_to_regcompass",
    "rc_compass_two_step_lp",
    "rc_hard_min_capacity",
    "rc_q95_calibrate",
    "rc_q95_shrink",
    "rc_q95_bootstrap_diagnostics",
    ".rc_weighted_q95_calibrate",
    ".rc_weighted_gene_score",
    ".rc_equal_sample_weights",
    ".rc_run_regcompass_stratum",
    ".rc_merge_stratum_layer1",
    ".rc_pando_gene_confidence",
    ".rc_pando_reaction_confidence",
    ".rc_pando_reaction_confidence_matrix",
    "rc_layer2_confidence_matrix",
    "rc_align_layer2_evidence",
    ".rc_apply_used_metacell_ids",
    ".rc_layer2_penalty_engine",
    "rc_layer2_penalty",
    "rc_layer2_support_penalties",
    "rc_layer2_reaction_type",
    "rc_layer2_support_penalty_for_type",
    "rc_layer2_has_gpr",
    ".rc_meta_module_one_hop"
  )
  expect_false(any(vapply(retired, exists, logical(1), inherits = TRUE)))
})

test_that("meta-module expansion exposes no one-hop controls", {
  retired_arguments <- c("include_one_hop", "one_hop_max_metabolite_degree")
  expect_false(any(retired_arguments %in%
                     names(formals(rc_expand_meta_module_reactions))))
  expect_false(any(retired_arguments %in%
                     names(formals(.rc_expand_meta_module_reactions_core))))
})

test_that("Seurat stack uses valid minimum bounds and exact stack metadata", {
  description <- utils::packageDescription("RegCompassR")
  imports <- description$Imports %||% ""
  expect_match(imports, "SeuratObject (>= 4.1.4)", fixed = TRUE)
  expect_match(imports, "Seurat (>= 4.4.0)", fixed = TRUE)
  expect_match(imports, "Signac (>= 1.11.0)", fixed = TRUE)

  expected <- c(
    SeuratObject = description[[
      "Config/RegCompassR/SeuratObjectVersion"
    ]],
    Seurat = description[["Config/RegCompassR/SeuratVersion"]],
    Signac = description[["Config/RegCompassR/SignacVersion"]]
  )
  expect_identical(
    unname(expected),
    c("4.1.4", "4.4.0", "1.11.0")
  )
  observed <- vapply(
    names(expected),
    function(package) as.character(utils::packageVersion(package)),
    character(1)
  )
  expect_identical(observed, expected)
})
