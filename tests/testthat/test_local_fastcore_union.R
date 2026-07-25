test_that("meta-module merging remains a catalogue, not a GEM", {
  biological <- data.frame(
    group_id = c("C1|T", "C2|T"),
    sample_id = c("C1|T", "C2|T"),
    module_id = c("C1|T::GRN0001", "C2|T::GRN0001"),
    reaction_id = c("Rcore", "Rcontext"),
    is_core = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  artifact <- list(
    group_id = "condition_pooled",
    grn_meta_modules = list(
      sample_status = data.frame(),
      tf_peak_gene_all = data.frame(),
      tf_peak_gene_significant = data.frame(),
      metabolic_gene_nodes = data.frame(),
      metabolic_gene_edges = data.frame(),
      core_gene_reaction = biological[biological$is_core, , drop = FALSE],
      reaction_membership = biological,
      meta_module_summary = data.frame()
    )
  )

  merged <- .rc_merge_stratum_meta_modules(list(artifact))

  expect_setequal(
    merged$merged_reaction_membership$reaction_id,
    c("Rcore", "Rcontext")
  )
  expect_setequal(
    merged$merged_core_reactions$reaction_id,
    "Rcore"
  )
  expect_false(merged$is_gem)
  expect_false(merged$fastcore_applied)
  expect_identical(
    merged$merge_source,
    "deduplicated_biological_meta_module_reactions"
  )
  expect_false(any(grepl("union", names(merged), ignore.case = TRUE)))
  expect_false(any(grepl("fastcore", merged$merged_reaction_membership$inclusion_stage,
                         ignore.case = TRUE)))
})

test_that("Stage 3 does not invoke local FASTCORE", {
  body_text <- paste(deparse(body(.rc_build_condition_meta_modules)),
                     collapse = "\n")
  expect_false(grepl("complete_stratum_meta_modules", body_text, fixed = TRUE))
  expect_false(grepl("local_fastcore", body_text, fixed = TRUE))
  expect_true(grepl("none_at_meta_module_stage", body_text, fixed = TRUE))
})

test_that("only the medium-specific cache creates union GEMs", {
  cache_body <- paste(
    deparse(body(.rc_build_medium_specific_union_gem_cache)),
    collapse = "\n"
  )
  expect_true(grepl("medium_scenario", cache_body, fixed = TRUE))
  expect_true(grepl("MEDIUM_UNION_GEM", cache_body, fixed = TRUE))
  expect_true(grepl("medium_specific_union_gem", cache_body, fixed = TRUE))
  expect_true(grepl("rc_build_meta_module_gem", cache_body, fixed = TRUE))
})

test_that("obsolete local FASTCORE controls are not public API", {
  exports <- getNamespaceExports("RegCompassR")
  expect_false(any(grepl("local.*fastcore|fastcore.*local", exports,
                         ignore.case = TRUE)))
  expect_false("rc_build_meta_module_gem" %in% exports)
})
