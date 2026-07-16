rc_balanced_solver_available <- function() {
  requireNamespace("highs", quietly = TRUE) ||
    requireNamespace("Rglpk", quietly = TRUE) ||
    requireNamespace("gurobi", quietly = TRUE)
}

rc_balanced_test_solver <- function() {
  if (requireNamespace("highs", quietly = TRUE)) return("highs")
  if (requireNamespace("Rglpk", quietly = TRUE)) return("glpk")
  if (requireNamespace("gurobi", quietly = TRUE)) return("gurobi")
  "highs"
}

rc_balanced_forward_toy <- function() {
  S <- matrix(
    c(
       1, -1, -1,  0,
       0,  1,  1, -1
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("A_c", "B_c"),
      c("EX_A", "Rcore", "Ralt", "EX_B")
    )
  )
  rc_make_gem(
    S,
    lb = c(EX_A = 0, Rcore = 0, Ralt = 0, EX_B = 0),
    ub = c(EX_A = 1000, Rcore = 1000, Ralt = 1000, EX_B = 1000),
    reaction_meta = data.frame(
      reaction_id = colnames(S),
      role = c("exchange", "internal", "internal", "exchange"),
      role_source = "curated",
      stringsAsFactors = FALSE
    )
  )
}

test_that("equal-sample weights give every biological sample equal total mass", {
  sample_ids <- c("S1", "S1", "S1", "S2", "S2", "S3")
  weights <- .rc_equal_sample_weights(sample_ids)
  totals <- tapply(weights, sample_ids, sum)
  expect_equal(unname(totals), rep(1 / 3, 3), tolerance = 1e-12)
  expect_equal(sum(weights), 1, tolerance = 1e-12)
})

test_that("weighted gene score and Q95 retain matrix dimensions and finite bounds", {
  expression <- matrix(
    c(0, 0, 10, 2, 2, 2),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), c("S1_M1", "S1_M2", "S2_M1"))
  )
  weights <- .rc_equal_sample_weights(c("S1", "S1", "S2"))
  score <- .rc_weighted_gene_score(expression, weights)
  expect_identical(dim(score), dim(expression))
  expect_identical(dimnames(score), dimnames(expression))
  expect_true(all(score >= 0 & score <= 1))

  capacity <- rbind(R1 = c(0.1, 0.2, 1.0), R2 = c(1, 1, 1))
  calibrated <- .rc_weighted_q95_calibrate(capacity, weights)
  expect_identical(dim(calibrated$C_rel), dim(capacity))
  expect_true(all(calibrated$C_rel <= 1, na.rm = TRUE))
  expect_true(all(calibrated$Q$sample_balanced))
})

test_that("local FASTCORE completes each GRN meta-module before union", {
  skip_if_not(rc_balanced_solver_available())
  gem <- rc_balanced_forward_toy()
  membership <- data.frame(
    group_id = "C1|S1|T",
    condition = "C1",
    sample_id = "S1",
    cell_type = "T",
    module_id = "C1|S1|T::GRN0001",
    reaction_id = c("Rcore", "Ralt"),
    is_core = c(TRUE, FALSE),
    inclusion_stage = c("core", "subsystem"),
    stringsAsFactors = FALSE
  )
  core <- membership[membership$is_core, , drop = FALSE]
  completed <- .rc_complete_stratum_meta_modules(
    list(reaction_membership = membership, core_gene_reaction = core),
    gem,
    outdir = tempfile("local_fastcore_test_"),
    local_fastcore_args = list(
      solver = rc_balanced_test_solver(),
      save_models = FALSE,
      strict = TRUE
    )
  )
  support <- completed$completed_reaction_membership$reaction_id[
    completed$completed_reaction_membership$local_fastcore_support
  ]
  expect_setequal(support, c("EX_A", "EX_B"))
  expect_true(all(c("Rcore", "Ralt") %in%
                    completed$completed_reaction_membership$reaction_id))
  expect_equal(completed$summary$n_local_fastcore_support, 2)
})

test_that("global union deduplicates locally completed support reactions", {
  biological <- data.frame(
    group_id = "C1|S1|T",
    sample_id = "S1",
    module_id = "C1|S1|T::GRN0001",
    reaction_id = "Rcore",
    is_core = TRUE,
    stringsAsFactors = FALSE
  )
  completed <- rbind(
    transform(
      biological,
      biological_meta_module_member = TRUE,
      local_fastcore_support = FALSE,
      inclusion_stage = "core"
    ),
    data.frame(
      group_id = rep("C1|S1|T", 2),
      sample_id = rep("S1", 2),
      module_id = rep("C1|S1|T::GRN0001", 2),
      reaction_id = c("EX_A", "EX_B"),
      is_core = FALSE,
      biological_meta_module_member = FALSE,
      local_fastcore_support = TRUE,
      inclusion_stage = "local_fastcore_support",
      stringsAsFactors = FALSE
    )
  )
  artifact <- list(
    group_id = "C1|S1|T",
    grn_meta_modules = list(
      sample_status = data.frame(),
      tf_peak_gene_all = data.frame(),
      tf_peak_gene_significant = data.frame(),
      metabolic_gene_nodes = data.frame(),
      metabolic_gene_edges = data.frame(),
      core_gene_reaction = biological,
      reaction_membership = biological,
      meta_module_summary = data.frame(),
      local_completed_reaction_membership = completed,
      local_fastcore_summary = data.frame(),
      local_fastcore_diagnostics = data.frame(),
      local_fastcore_completion_iterations = data.frame()
    )
  )
  merged <- .rc_merge_stratum_meta_modules(list(artifact))
  expect_setequal(
    merged$global_reaction_membership$reaction_id,
    c("Rcore", "EX_A", "EX_B")
  )
  expect_equal(
    merged$global_union_source,
    "deduplicated_local_fastcore_completed_meta_modules"
  )
  expect_true(all(
    merged$global_reaction_membership$inclusion_stage[
      merged$global_reaction_membership$reaction_id %in% c("EX_A", "EX_B")
    ] == "global_union_local_fastcore_support"
  ))
})

test_that("limma correction refuses to remove biological sample IDs", {
  skip_if_not_installed("limma")
  expression <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("M1", "M2"))
  )
  meta <- data.frame(
    sample_id = c("S1", "S2"),
    condition = c("A", "B"),
    cell_type = c("T", "T"),
    stringsAsFactors = FALSE
  )
  params <- .rc_normalize_calibration_params(
    list(
      expression_batch_correction = "limma",
      technical_batch_cols = "sample_id",
      preserve_design_cols = "condition"
    ),
    sample_col = "sample_id",
    condition_col = "condition",
    celltype_col = "cell_type"
  )
  expect_error(
    .rc_apply_limma_batch_correction(
      expression,
      meta,
      params,
      sample_col = "sample_id"
    ),
    "biological replicate"
  )
})
