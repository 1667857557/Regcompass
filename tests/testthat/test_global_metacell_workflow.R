test_that("integrated workflow keeps one barrier between the two worker pools", {
  core_text <- paste(deparse(body(.rc_run_regcompass_engine)), collapse = "\n")
  wrapper_text <- paste(deparse(body(rc_run_regcompass)), collapse = "\n")
  expect_match(core_text, ".rc_run_regcompass_stratum", fixed = TRUE)
  expect_match(core_text, "upstream_complete_barrier", fixed = TRUE)
  expect_match(core_text, ".rc_release_bpparam(upstream_param)", fixed = TRUE)
  expect_match(core_text, ".rc_merge_stratum_meta_modules", fixed = TRUE)
  expect_match(core_text, ".rc_merge_stratum_layer1", fixed = TRUE)
  expect_match(wrapper_text, "inference_unit", fixed = TRUE)
  expect_match(wrapper_text, "04_model_cache", fixed = TRUE)
  expect_equal(
    eval(formals(rc_run_regcompass)$inference_unit),
    c("sample_celltype", "metacell")
  )
})


test_that("integrated workflow uses fixed metacell gamma and skips low-yield strata", {
  workflow_text <- paste(deparse(body(.rc_run_regcompass_engine)), collapse = "\n")
  worker_text <- paste(
    deparse(body(.rc_run_regcompass_stratum)),
    collapse = "\n"
  )
  metacell_formals <- names(formals(rc_make_supercell2_metacells))
  expect_false("adaptive_gamma" %in% metacell_formals)
  expect_false(grepl("adaptive_gamma", workflow_text, fixed = TRUE))
  expect_false(grepl("adaptive_gamma", worker_text, fixed = TRUE))
  expect_match(worker_text, "skipped_too_few_metacells", fixed = TRUE)
  expect_match(workflow_text, "n_skipped_too_few_metacells", fixed = TRUE)
  expect_match(worker_text, "or_method", fixed = TRUE)
})


test_that("strict-stratum worker runs Pando and local FASTCORE", {
  body_text <- paste(
    deparse(body(.rc_run_regcompass_stratum)),
    collapse = "\n"
  )
  expect_match(body_text, "rc_make_supercell2_metacells", fixed = TRUE)
  expect_match(body_text, "rc_run_pando_meta_modules", fixed = TRUE)
  expect_match(body_text, ".rc_complete_stratum_meta_modules", fixed = TRUE)
  expect_false(grepl("LinkPeaks", body_text, fixed = TRUE))
})


test_that("global Layer 1 uses equal total sample weight", {
  artifacts <- list(
    list(
      capacity_params = list(
        promiscuity_mode = "sqrt", and_method = "boltzmann",
        tau = 0.20, or_method = "sum_sqrtK"
      ),
      calibration_params = list(sample_balance = TRUE),
      layer1 = list(
        rna_metacell_logcpm = matrix(1, 1, 2, dimnames = list("G1", c("u1", "u2"))),
        reaction_confidence = matrix(0.8, 1, 2, dimnames = list("R1", c("u1", "u2"))),
        unit_meta = data.frame(
          pool_id = c("u1", "u2"), unit_id = c("u1", "u2"),
          sample_id = c("S1", "S1"), condition = "A", cell_type = "T"
        )
      )
    ),
    list(
      capacity_params = list(
        promiscuity_mode = "sqrt", and_method = "boltzmann",
        tau = 0.20, or_method = "sum_sqrtK"
      ),
      calibration_params = list(sample_balance = TRUE),
      layer1 = list(
        rna_metacell_logcpm = matrix(3, 1, 1, dimnames = list("G1", "u3")),
        reaction_confidence = matrix(0.4, 1, 1, dimnames = list("R1", "u3")),
        unit_meta = data.frame(
          pool_id = "u3", unit_id = "u3", sample_id = "S2",
          condition = "B", cell_type = "T"
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
  totals <- tapply(
    as.numeric(out$sample_balance_weights),
    out$unit_meta$sample_id,
    sum
  )
  expect_equal(as.numeric(totals), c(0.5, 0.5), tolerance = 1e-12)
  expect_match(out$capacity_calibration_scope, "equal_sample_weighted", fixed = TRUE)
})


test_that("global union contains biological and local support reactions", {
  artifact <- list(
    group_id = "A|S1|T",
    grn_meta_modules = list(
      sample_status = data.frame(status = "ok"),
      tf_peak_gene_all = data.frame(),
      tf_peak_gene_significant = data.frame(),
      metabolic_gene_nodes = data.frame(),
      metabolic_gene_edges = data.frame(),
      core_gene_reaction = data.frame(
        sample_id = "S1", module_id = "M1", reaction_id = "R1", is_core = TRUE
      ),
      reaction_membership = data.frame(
        sample_id = "S1", module_id = "M1", reaction_id = c("R1", "R2")
      ),
      local_completed_reaction_membership = data.frame(
        sample_id = "S1", module_id = "M1", reaction_id = c("R1", "R2", "R3")
      ),
      meta_module_summary = data.frame()
    )
  )
  out <- .rc_merge_stratum_meta_modules(list(artifact))
  expect_setequal(out$global_reaction_membership$reaction_id, c("R1", "R2", "R3"))
  expect_equal(
    out$global_reaction_membership$inclusion_stage[
      out$global_reaction_membership$reaction_id == "R3"
    ],
    "global_union_local_fastcore_support"
  )
})
