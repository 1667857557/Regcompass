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

test_that("structural override flags retain reaction identifiers", {
  E <- matrix(
    0.1,
    nrow = 2,
    dimnames = list(c("EX_A", "R_B"), "mc1")
  )
  roles <- data.frame(
    reaction_id = rownames(E),
    role = c("exchange", "internal"),
    role_source = c("id_pattern", "curated"),
    stringsAsFactors = FALSE
  )
  answer <- rc_compute_multiome_penalty(E, reaction_roles = roles)
  expect_identical(names(answer$components$role_override_flag), rownames(E))
  expect_true(answer$components$role_override_flag[["EX_A"]])
  expect_false(answer$components$role_override_flag[["R_B"]])
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

test_that("condition pooling requires replicated conditions in strict mode", {
  meta <- data.frame(
    sample_id = c("A1", "A2", "B1", "B2"),
    condition = c("A", "A", "B", "B"),
    cell_type = "T",
    stringsAsFactors = FALSE
  )
  design <- .rc_condition_pool_design_summary(
    meta,
    "sample_id",
    "condition",
    "cell_type",
    strict_biological_defaults = TRUE
  )
  expect_equal(design$condition_sample_count$n_biological_samples, c(2, 2))

  one_sample <- meta[meta$sample_id != "B2", , drop = FALSE]
  expect_error(
    .rc_condition_pool_design_summary(
      one_sample,
      "sample_id",
      "condition",
      "cell_type",
      strict_biological_defaults = TRUE
    ),
    "at least two biological samples"
  )
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

test_that("meta-module adds one bounded non-structural metabolite hop", {
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
  expect_setequal(expanded$reaction_membership$reaction_id, c("R1", "R2"))
  stage <- stats::setNames(
    expanded$reaction_membership$inclusion_stage,
    expanded$reaction_membership$reaction_id
  )
  expect_identical(stage[["R2"]], "one_hop_metabolite_neighbor")
  expect_equal(expanded$summary$n_one_hop_added, 1)
  expect_false(any(c("R3", "EX_M1") %in%
                     expanded$reaction_membership$reaction_id))
})
