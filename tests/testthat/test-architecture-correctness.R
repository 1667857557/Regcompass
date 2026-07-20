test_that("absolute RNA support preserves zero and constant abundance", {
  zero <- matrix(0, nrow = 1, ncol = 3,
                 dimnames = list("g0", paste0("m", 1:3)))
  high <- matrix(log1p(10), nrow = 1, ncol = 3,
                 dimnames = list("g1", paste0("m", 1:3)))

  expect_equal(as.numeric(rc_gene_score(zero)), rep(0, 3))
  expected_high <- log1p(10) / (log1p(10) + 1)
  expect_equal(as.numeric(rc_gene_score(high)), rep(expected_high, 3))
})

test_that("nested GPR rules preserve Boolean logic", {
  parsed <- rc_parse_gpr_simple("g1 and (g2 or g3)")
  keys <- sort(vapply(parsed, function(x) paste(sort(x), collapse = "+"),
                      character(1)))
  expect_equal(keys, c("g1+g2", "g1+g3"))

  parsed2 <- rc_parse_gpr_simple("(g1 and g2) or (g3 and (g4 or g5))")
  keys2 <- sort(vapply(parsed2, function(x) paste(sort(x), collapse = "+"),
                       character(1)))
  expect_equal(keys2, c("g1+g2", "g3+g4", "g3+g5"))
  expect_error(rc_parse_gpr_simple("g1 and (g2 or)"), "Malformed GPR")
})

test_that("TF-ATAC integration is zero preserving and signed", {
  C <- matrix(c(0, 0.2, 0.2), nrow = 3,
              dimnames = list(c("zero", "activated", "repressed"), "u1"))
  R <- matrix(c(1, 1, -1), nrow = 3, dimnames = dimnames(C))
  out <- .rc_integrate_regulatory_support_v170(C, R, alpha = 1)

  expect_equal(out["zero", "u1"], 0)
  expect_gt(out["activated", "u1"], C["activated", "u1"])
  expect_lt(out["repressed", "u1"], C["repressed", "u1"])
  expect_true(all(out >= 0 & out <= 1))
})

test_that("COMPASS-like penalty is positive and monotonically decreasing", {
  expression <- matrix(
    c(0, 1, 3, NA_real_),
    nrow = 4,
    dimnames = list(c("zero", "low", "high", "missing"), "u1")
  )
  answer <- rc_compute_multiome_penalty(expression)
  penalty <- answer$penalty[, "u1"]

  expect_equal(penalty[["zero"]], 1)
  expect_equal(penalty[["missing"]], 1)
  expect_gt(penalty[["zero"]], penalty[["low"]])
  expect_gt(penalty[["low"]], penalty[["high"]])
  expect_true(all(is.finite(penalty) & penalty > 0))
  expect_identical(answer$evidence_policy, "penalty_only")
  expect_identical(answer$penalty_version,
                   "v1.7.0_gene_integrated_multiome_penalty")
})

test_that("shared-TF projection retains regulator and sign metadata", {
  edges <- data.frame(
    sample_id = c("s1", "s1"),
    tf = c("TF1", "TF1"),
    target = c("G1", "G2"),
    estimate = c(1, -2),
    stringsAsFactors = FALSE
  )
  projected <- rc_project_metabolic_grn(
    edges, metabolic_genes = c("G1", "G2"),
    top_k = 5, min_shared_tfs = 1, min_tf_jaccard = 0
  )
  expect_equal(nrow(projected$edges), 1)
  expect_equal(projected$edges$regulator_set, "TF1")
  expect_equal(projected$edges$regulatory_relation, "discordant")
  expect_lt(projected$edges$signed_projection_weight, 0)
  expect_true(projected$edges$direction_and_sign_preserved)
})

test_that("condition-pooled metacell is the canonical inference unit", {
  expect_identical(eval(formals(rc_run_regcompass)$inference_unit), "metacell")
  expect_identical(eval(formals(rc_run_regcompass_one_shot)$medium_scenario),
                   "physiologic")
})
