test_that("integrated workflow fuses strict-stratum upstream work", {
  body_text <- paste(deparse(body(rc_run_regcompass)), collapse = "\n")
  expect_match(body_text, ".rc_run_regcompass_stratum", fixed = TRUE)
  expect_match(body_text, ".rc_merge_stratum_meta_modules", fixed = TRUE)
  expect_match(body_text, ".rc_release_bpparam(upstream_param)", fixed = TRUE)
  expect_match(body_text, "unit = \"metacell\"", fixed = TRUE)
  expect_match(body_text, "global_reaction_membership", fixed = TRUE)
})

test_that("global meta-module union preserves source tables and creates one canonical module", {
  artifact <- list(
    group_id = "A|S1|T",
    grn_meta_modules = list(
      sample_status = data.frame(group_id = "A|S1|T", status = "ok"),
      tf_peak_gene_all = data.frame(),
      tf_peak_gene_significant = data.frame(),
      metabolic_gene_nodes = data.frame(),
      metabolic_gene_edges = data.frame(),
      core_gene_reaction = data.frame(
        sample_id = "S1", module_id = "A|S1|T::GRN0001",
        reaction_id = c("R1", "R2"), is_core = c(TRUE, FALSE),
        stringsAsFactors = FALSE
      ),
      reaction_membership = data.frame(
        sample_id = "S1", module_id = "A|S1|T::GRN0001",
        reaction_id = c("R1", "R2", "R3"),
        stringsAsFactors = FALSE
      ),
      meta_module_summary = data.frame()
    )
  )
  out <- .rc_merge_stratum_meta_modules(list(artifact))
  expect_equal(out$global_core_reactions$reaction_id, "R1")
  expect_setequal(out$global_reaction_membership$reaction_id, c("R1", "R2", "R3"))
  expect_true(all(out$global_reaction_membership$sample_id == "global"))
  expect_true(all(out$global_reaction_membership$module_id == "GLOBAL_UNION"))
})

test_that("microCOMPASS applies every shared model target to every unit", {
  body_text <- paste(deparse(body(rc_run_microcompass)), collapse = "\n")
  expect_match(body_text, "expand.grid(row_id = row_ids, unit_id = units", fixed = TRUE)
  expect_false(grepl("unit_sample ==", body_text, fixed = TRUE))
  expect_match(body_text, "shared_gem_scope = \"all_metacells\"", fixed = TRUE)
})

test_that("legacy integrated and sample-specific cache APIs are retired", {
  exports <- getNamespaceExports("RegCompassR")
  expect_false("rc_run_regcompass_multiome_metacell" %in% exports)
  expect_false("rc_build_meta_module_gem_cache" %in% exports)
  expect_false("rc_load_metacell_object_from_run" %in% exports)
})
