test_that("v1.7.0 regulatory integration is bounded and zero preserving", {
  C <- matrix(
    c(0, 0.20, 0.20, 0.80),
    nrow = 4,
    dimnames = list(paste0("g", 1:4), "mc1")
  )
  R <- matrix(
    c(1, 1, -1, 1),
    nrow = 4,
    dimnames = dimnames(C)
  )
  out <- .rc_integrate_regulatory_support_v170(C, R, alpha = 1)

  expect_equal(out["g1", "mc1"], 0)
  expect_gt(out["g2", "mc1"], C["g2", "mc1"])
  expect_lt(out["g3", "mc1"], C["g3", "mc1"])
  expect_gt(out["g4", "mc1"], C["g4", "mc1"])
  expect_true(all(out >= 0 & out <= 1))
})

test_that("v1.7.0 reaction penalty is positive and decreases with expression", {
  E <- matrix(
    c(0, 1, 3, NA_real_),
    nrow = 4,
    dimnames = list(paste0("R", 1:4), "mc1")
  )
  answer <- rc_compute_multiome_penalty(E)
  P <- answer$penalty[, "mc1"]

  expect_equal(P[["R1"]], 1)
  expect_gt(P[["R1"]], P[["R2"]])
  expect_gt(P[["R2"]], P[["R3"]])
  expect_equal(P[["R4"]], 1)
  expect_true(all(is.finite(P) & P > 0))
  expect_identical(answer$penalty_formula, "1 / (1 + log2(1 + E_multiome))")
  expect_identical(answer$evidence_policy, "penalty_only")
  expect_false("P_conf" %in% names(answer$components))
})

test_that("only structural reactions receive fixed penalty overrides", {
  E <- matrix(
    c(0.1, 0.1, 0.1),
    nrow = 3,
    dimnames = list(c("EX_A", "T_A", "R_B"), "mc1")
  )
  roles <- data.frame(
    reaction_id = rownames(E),
    role = c("exchange", "transport", "internal"),
    role_source = c("id_pattern", "curated", "curated"),
    stringsAsFactors = FALSE
  )
  answer <- rc_compute_multiome_penalty(E, reaction_roles = roles)
  flag <- answer$components$role_override_flag
  expect_identical(names(flag), rownames(E))
  expect_true(flag[["EX_A"]])
  expect_false(flag[["T_A"]])
  expect_false(flag[["R_B"]])
  expect_equal(
    answer$penalty["T_A", "mc1"],
    1 / (1 + log2(1 + E["T_A", "mc1"]))
  )
})

test_that("condition-pooled grouping excludes biological sample", {
  expect_identical(
    .rc_condition_group_cols("condition", "cell_type"),
    c("condition", "cell_type")
  )
  expect_error(
    .rc_condition_group_cols("condition", "condition"),
    "must be distinct"
  )
})

test_that("one biological sample per condition is accepted descriptively", {
  meta <- data.frame(
    sample_id = c("A1", "B1"),
    condition = c("A", "B"),
    cell_type = "T",
    stringsAsFactors = FALSE
  )
  expect_warning(
    design <- .rc_condition_pool_design_summary(
      meta,
      "sample_id",
      "condition",
      "cell_type",
      strict_biological_defaults = FALSE
    ),
    "at least two biological samples"
  )
  expect_equal(design$condition_sample_count$n_biological_samples, c(1, 1))
  expect_equal(nrow(design$low_replication_strata), 2L)
  expect_false("strict_biological_defaults" %in% names(formals(rc_run_regcompass)))
  expect_false("inference_unit" %in% names(formals(rc_run_regcompass)))
})

test_that("per-metacell regulatory state uses ATAC rather than TF RNA", {
  body_text <- paste(
    deparse(body(.rc_condition_gene_regulatory_modifier)),
    collapse = "\n"
  )
  expect_match(body_text, ".rc_pando_assay_data(object, atac_assay)", fixed = TRUE)
  expect_false(grepl("tf_score", body_text, fixed = TRUE))
  expect_false(grepl("rna_assay", body_text, fixed = TRUE))
})

test_that("single-condition scoring uses penalty per required target flux", {
  row_ids <- c(
    "reaction=R1::direction=forward::medium=base",
    "reaction=R2::direction=forward::medium=base"
  )
  units <- c("u1", "u2")
  microcompass <- list(
    penalty = matrix(
      c(0.2, 0.4, 0.3, 0.5),
      nrow = 2,
      dimnames = list(row_ids, units)
    ),
    vmax = matrix(
      c(10, 1, 10, 1),
      nrow = 2,
      dimnames = list(row_ids, units)
    ),
    unit_meta = data.frame(
      unit_id = units,
      condition = "A",
      cell_type = "T",
      stringsAsFactors = FALSE
    ),
    params = list(omega = 0.95)
  )
  answer <- .rc_condition_penalty_comparison(microcompass)
  expect_identical(answer$analysis_mode, "single_condition_reaction_ranking")
  expect_identical(answer$ranking_formula, "penalty / (omega * vmax)")
  expect_equal(nrow(answer$ranking), 2L)
  expect_equal(answer$ranking$reaction_id, c("R1", "R2"))
  expect_equal(answer$ranking$priority_rank, c(1L, 2L))
  expect_equal(
    answer$ranking$median_penalty_per_target_flux,
    c(0.25 / 9.5, 0.45 / 0.95)
  )
  expect_equal(nrow(answer$contrast), 0L)
})

test_that("multiple conditions produce every pairwise descriptive comparison", {
  row_id <- "reaction=R1::direction=forward::medium=base"
  units <- c("uA", "uB", "uC")
  microcompass <- list(
    penalty = matrix(
      c(0.5, 0.3, 0.2),
      nrow = 1,
      dimnames = list(row_id, units)
    ),
    vmax = matrix(
      2,
      nrow = 1,
      ncol = 3,
      dimnames = list(row_id, units)
    ),
    unit_meta = data.frame(
      unit_id = units,
      condition = c("A", "B", "C"),
      cell_type = "T",
      stringsAsFactors = FALSE
    ),
    params = list(omega = 0.95)
  )
  answer <- .rc_condition_penalty_comparison(microcompass)
  expect_identical(
    answer$analysis_mode,
    "multi_condition_reaction_ranking_and_pairwise_comparison"
  )
  expect_equal(nrow(answer$ranking), 3L)
  expect_equal(nrow(answer$contrast), 3L)
  expect_setequal(
    paste(answer$contrast$condition_a, answer$contrast$condition_b, sep = "-"),
    c("A-B", "A-C", "B-C")
  )
})

test_that("shared-model ranking rejects unit-dependent vmax", {
  row_id <- "reaction=R1::direction=forward::medium=base"
  microcompass <- list(
    penalty = matrix(
      c(0.5, 0.3), nrow = 1,
      dimnames = list(row_id, c("u1", "u2"))
    ),
    vmax = matrix(
      c(1, 2), nrow = 1,
      dimnames = list(row_id, c("u1", "u2"))
    ),
    unit_meta = data.frame(
      unit_id = c("u1", "u2"),
      condition = c("A", "B"),
      cell_type = "T",
      stringsAsFactors = FALSE
    ),
    params = list(omega = 0.95)
  )
  expect_error(
    .rc_condition_penalty_comparison(microcompass),
    "vmax differs across metacells"
  )
})

test_that("meta-module expansion excludes metabolite-neighbour reactions", {
  S <- matrix(
    c(
      -1,  0,  0,  1,
       1, -1,  0,  0,
       0,  1, -1,  0,
       0,  0,  1,  0
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(
      paste0("M", 1:4),
      c("R1", "R2", "R3", "EX_M1")
    )
  )
  reaction_meta <- data.frame(
    reaction_id = colnames(S),
    subsystem = c("A", "B", "C", "D"),
    metabolic_module = c("A", "B", "C", "D"),
    role = c("internal", "internal", "internal", "exchange"),
    role_source = "curated",
    stringsAsFactors = FALSE
  )
  gem <- rc_make_gem(
    S,
    lb = rep(0, ncol(S)),
    ub = rep(1000, ncol(S)),
    reaction_meta = reaction_meta
  )
  core <- data.frame(
    sample_id = "A",
    module_id = "A::GRN0001",
    gene = "G1",
    reaction_id = "R1",
    stringsAsFactors = FALSE
  )
  expanded <- rc_expand_meta_module_reactions(gem, core)

  expect_identical(
    unique(as.character(expanded$reaction_membership$reaction_id)),
    "R1"
  )
  expect_false(any(grepl("one_hop", names(expanded$summary), fixed = TRUE)))
  expect_false(any(c(
    "include_one_hop", "one_hop_max_metabolite_degree"
  ) %in% names(formals(rc_expand_meta_module_reactions))))
  expect_false(any(c(
    "include_one_hop", "one_hop_max_metabolite_degree"
  ) %in% names(formals(.rc_expand_meta_module_reactions_core))))
  expect_false(exists(".rc_meta_module_one_hop", inherits = TRUE))
})
