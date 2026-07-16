test_that("integrated workflow has an all-strata barrier before global work", {
  body_text <- paste(
    deparse(body(rc_run_regcompass_balanced_v15)),
    collapse = "\n"
  )
  expect_match(body_text, ".rc_run_regcompass_stratum", fixed = TRUE)
  expect_match(body_text, "upstream_complete_barrier", fixed = TRUE)
  expect_match(body_text, "not all retained strata completed successfully", fixed = TRUE)
  expect_match(body_text, ".rc_release_bpparam(upstream_param)", fixed = TRUE)
  expect_match(body_text, ".rc_merge_stratum_meta_modules", fixed = TRUE)
  expect_match(body_text, ".rc_merge_stratum_layer1", fixed = TRUE)
  expect_match(body_text, "unit = \"metacell\"", fixed = TRUE)
  expect_match(body_text, "global_reaction_membership", fixed = TRUE)

  wrapper_text <- paste(deparse(body(rc_run_regcompass)), collapse = "\n")
  expect_match(wrapper_text, "sample_balanced_gene_score", fixed = TRUE)
  expect_match(wrapper_text, "local_fastcore_before_global_union", fixed = TRUE)
})

test_that("global Layer 1 recomputes sample-balanced gene scores and uses Pando confidence", {
  artifacts <- list(
    list(
      capacity_params = list(
        promiscuity_mode = "sqrt", and_method = "boltzmann",
        tau = 0.20, or_method = "sum_sqrtK"
      ),
      calibration_params = list(
        sample_balance = TRUE,
        sample_balance_col = "sample_id",
        expression_batch_correction = "none",
        technical_batch_cols = character(),
        preserve_design_cols = c("condition", "cell_type")
      ),
      layer1 = list(
        rna_metacell_logcpm = matrix(
          1, nrow = 1, dimnames = list("G1", "u1")
        ),
        reaction_confidence = matrix(
          0.8, nrow = 1, dimnames = list("R1", "u1")
        ),
        unit_meta = data.frame(
          pool_id = "u1", unit_id = "u1", sample_id = "S1",
          condition = "A", cell_type = "T", stringsAsFactors = FALSE
        )
      )
    ),
    list(
      capacity_params = list(
        promiscuity_mode = "sqrt", and_method = "boltzmann",
        tau = 0.20, or_method = "sum_sqrtK"
      ),
      calibration_params = list(
        sample_balance = TRUE,
        sample_balance_col = "sample_id",
        expression_batch_correction = "none",
        technical_batch_cols = character(),
        preserve_design_cols = c("condition", "cell_type")
      ),
      layer1 = list(
        rna_metacell_logcpm = matrix(
          3, nrow = 1, dimnames = list("G1", "u2")
        ),
        reaction_confidence = matrix(
          0.4, nrow = 1, dimnames = list("R1", "u2")
        ),
        unit_meta = data.frame(
          pool_id = "u2", unit_id = "u2", sample_id = "S2",
          condition = "B", cell_type = "T", stringsAsFactors = FALSE
        )
      )
    )
  )
  gem <- list(gpr_table = data.frame(
    reaction_id = "R1", and_group_id = 1, gene = "G1"
  ))
  out <- .rc_merge_stratum_layer1(
    artifacts, gem, single_cell_genes = "G1",
    sample_col = "sample_id", condition_col = "condition",
    celltype_col = "cell_type"
  )
  expect_equal(
    out$capacity_calibration_scope,
    "equal_sample_weighted_none_global_gene_score_and_reaction_q95"
  )
  expect_equal(
    out$reaction_confidence_source,
    "pando_internal_peak_gene_accessibility"
  )
  expect_equal(
    unique(as.character(out$q95_diagnostics$stratum)),
    "global_sample_balanced"
  )
  expect_true(out$calibration_params$sample_balance)
  expect_equal(
    as.numeric(tapply(
      out$sample_balance_weights,
      out$unit_meta$sample_id,
      sum
    )),
    c(0.5, 0.5)
  )
  expect_identical(colnames(out$C_rel), c("u1", "u2"))
  expect_equal(out$reaction_confidence["R1", ], c(u1 = 0.8, u2 = 0.4))
  expect_lt(out$global_gene_score["G1", "u1"], out$global_gene_score["G1", "u2"])
  expect_lt(out$C_raw["R1", "u1"], out$C_raw["R1", "u2"])
})

test_that("medium normalization preserves an existing no-constraint marker", {
  base <- .rc_normalize_medium_scenarios(NULL)
  expect_true(base$.no_constraints)
  normalized_again <- .rc_normalize_medium_scenarios(base)
  expect_true(normalized_again$.no_constraints)
})

test_that("integrated stratum worker uses Pando without rerunning LinkPeaks", {
  body_text <- paste(
    deparse(body(.rc_run_regcompass_stratum_v14)),
    collapse = "\n"
  )
  expect_match(body_text, "rc_run_pando_meta_modules", fixed = TRUE)
  expect_match(
    body_text,
    "pando_internal_peak_gene_accessibility",
    fixed = TRUE
  )
  expect_false(grepl("rc_run_layer1_from_metacells", body_text, fixed = TRUE))
  expect_false(grepl("Signac::LinkPeaks", body_text, fixed = TRUE))

  completion_text <- paste(
    deparse(body(.rc_run_regcompass_stratum_balanced_v15)),
    collapse = "\n"
  )
  expect_match(
    completion_text,
    ".rc_complete_stratum_meta_modules",
    fixed = TRUE
  )
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
      local_completed_reaction_membership = data.frame(
        sample_id = "S1", module_id = "A|S1|T::GRN0001",
        reaction_id = c("R1", "R2", "R3", "SUPPORT"),
        local_fastcore_support = c(FALSE, FALSE, FALSE, TRUE),
        stringsAsFactors = FALSE
      ),
      meta_module_summary = data.frame(),
      local_fastcore_summary = data.frame(),
      local_fastcore_diagnostics = data.frame(),
      local_fastcore_completion_iterations = data.frame()
    )
  )
  out <- .rc_merge_stratum_meta_modules(list(artifact))
  expect_equal(out$global_core_reactions$reaction_id, "R1")
  expect_setequal(
    out$global_reaction_membership$reaction_id,
    c("R1", "R2", "R3", "SUPPORT")
  )
  expect_true(all(out$global_reaction_membership$sample_id == "global"))
  expect_true(all(out$global_reaction_membership$module_id == "GLOBAL_UNION"))
  expect_equal(
    out$global_reaction_membership$inclusion_stage[
      out$global_reaction_membership$reaction_id == "SUPPORT"
    ],
    "global_union_local_fastcore_support"
  )
})

test_that("microCOMPASS applies each shared model to every metacell task", {
  body_text <- paste(deparse(body(rc_run_microcompass)), collapse = "\n")
  expect_match(body_text, "expand.grid(model_key = unique_model_keys, unit_id = units", fixed = TRUE)
  expect_match(body_text, "run_one_metacell", fixed = TRUE)
  expect_false(grepl("unit_sample ==", body_text, fixed = TRUE))
  expect_match(body_text, "shared_gem_scope = \"all_metacells\"", fixed = TRUE)
  expect_match(body_text, "parallel_task = \"shared_model_by_metacell\"", fixed = TRUE)
})

test_that("legacy integrated and sample-specific cache APIs are retired", {
  exports <- getNamespaceExports("RegCompassR")
  expect_false("rc_run_regcompass_multiome_metacell" %in% exports)
  expect_false("rc_build_meta_module_gem_cache" %in% exports)
  expect_false("rc_load_metacell_object_from_run" %in% exports)
})
