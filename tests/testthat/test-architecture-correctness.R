test_that("absolute gene evidence preserves zero and constant abundance", {
  zero <- matrix(0, nrow = 1, ncol = 3,
                 dimnames = list("g0", paste0("m", 1:3)))
  high <- matrix(log1p(10), nrow = 1, ncol = 3,
                 dimnames = list("g1", paste0("m", 1:3)))

  expect_equal(as.numeric(rc_gene_score(zero)), rep(0, 3))
  expected_high <- log1p(10) / (log1p(10) + 1)
  expect_equal(as.numeric(rc_gene_score(high)), rep(expected_high, 3))
  expect_true(all(rc_gene_score(high) > 0.5))

  raw <- rbind(zero = c(0, 0, 0), constant_high = c(0.8, 0.8, 0.8))
  colnames(raw) <- paste0("m", 1:3)
  calibrated <- rc_q95_calibrate(
    raw,
    bootstrap = FALSE,
    weights = rep(1 / 3, 3)
  )
  expect_equal(as.numeric(calibrated$C_rel["zero", ]), rep(0, 3))
  expect_equal(as.numeric(calibrated$C_rel["constant_high", ]), rep(0.8, 3))
  expect_true(all(is.na(calibrated$C_within_reaction_relative)))
  expect_true(all(calibrated$Q$calibration_role ==
                    "diagnostic_only_not_lp_capacity"))
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

test_that("canonical GPR defaults preserve bottlenecks and avoid count bias", {
  previous <- options(RegCompassR.strict_gpr_defaults = TRUE)
  on.exit(options(previous), add = TRUE)

  gene <- matrix(c(0, 1), nrow = 2, ncol = 1,
                 dimnames = list(c("g1", "g2"), "m1"))
  gpr <- list(
    complex = list(c("g1", "g2")),
    isoenzyme = list(c("g1"), c("g2")),
    extra_annotation = list(c("g2"))
  )
  capacity <- rc_reaction_capacity(gpr, gene)
  expect_equal(capacity["complex", "m1"], 0)
  expect_equal(capacity["isoenzyme", "m1"], 1)
  expect_equal(capacity["extra_annotation", "m1"], 1)
})

test_that("regulatory support is neutral when missing and penalizes repression", {
  capacity <- matrix(0.8, nrow = 1, ncol = 4,
                     dimnames = list("r1", paste0("u", 1:4)))
  regulation <- matrix(c(NA, 0.5, 0.25, 0.75), nrow = 1,
                       dimnames = dimnames(capacity))
  answer <- rc_compute_multiome_penalty(capacity, regulation)

  expect_equal(answer$penalty["r1", "u1"],
               answer$penalty["r1", "u2"], tolerance = 1e-8)
  expect_gt(answer$penalty["r1", "u3"],
            answer$penalty["r1", "u2"])
  expect_equal(answer$penalty["r1", "u4"],
               answer$penalty["r1", "u2"], tolerance = 1e-8)
  expect_true(answer$components$missing_regulatory_support_flag["r1", "u1"])
})

test_that("relative penalty rank is stable and explicitly not probability", {
  penalty <- rbind(variable = c(1, 1, 1.0001), constant = c(2, 2, 2))
  colnames(penalty) <- paste0("u", seq_len(ncol(penalty)))
  feasible <- matrix(TRUE, nrow = 2, ncol = 3,
                     dimnames = dimnames(penalty))
  score <- rc_compass_score_from_penalty(penalty, feasible)

  expect_true(all(score["variable", ] >= 0 & score["variable", ] <= 1))
  expect_true(all(is.na(score["constant", ])))
  expect_identical(
    attr(score, "score_semantics"),
    "within_target_relative_penalty_rank_not_probability"
  )
  expect_true(attr(score, "noninformative_target")[["constant"]])
})

test_that("Layer 2 feasibility is aligned by identifiers", {
  penalty <- matrix(
    c(1, 4, 2, 5, 3, 6), nrow = 2,
    dimnames = list(c("r1", "r2"), c("u1", "u2", "u3"))
  )
  feasible_aligned <- matrix(
    c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE), nrow = 2,
    dimnames = dimnames(penalty)
  )
  feasible <- feasible_aligned[c("r2", "r1"), c("u3", "u1", "u2")]
  score <- rc_compass_score_from_penalty(penalty, feasible)
  expect_true(is.na(score["r1", "u2"]))
  expect_true(all(is.finite(score["r1", c("u1", "u3")])))
  expect_error(
    rc_compass_score_from_penalty(penalty, feasible[, "u1", drop = FALSE]),
    "identical target and unit IDs",
    fixed = TRUE
  )
})

test_that("permissive medium is labelled as a technical baseline", {
  S <- Matrix::Matrix(matrix(c(-1, 1), nrow = 1), sparse = TRUE)
  rownames(S) <- "m_e"
  colnames(S) <- c("EX_m", "R1")
  gem <- list(
    S = S,
    lb = stats::setNames(c(-1000, 0), colnames(S)),
    ub = stats::setNames(c(1000, 1000), colnames(S)),
    reaction_meta = data.frame(
      reaction_id = colnames(S),
      role = c("exchange", "internal"),
      stringsAsFactors = FALSE
    )
  )

  medium <- rc_make_medium_scenarios(
    gem, scenario = "permissive_all_exchange"
  )
  expect_true(all(medium$medium_scenario_id == "permissive_all_exchange"))
  expect_true(all(medium$assumption_level == "technical_sensitivity_baseline"))
  expect_true(all(!medium$concentration_used_for_rate_bound))

  annotated <- gem
  annotated$reaction_meta$reaction_name <- c("Exchange of glucose", "internal")
  named_medium <- rc_make_medium_scenarios(
    annotated,
    scenario = "normal_human_plasma",
    strict_preset_matching = FALSE
  )
  expect_true(all(named_medium$medium_scenario_id == "normal_human_plasma"))
  expect_true(all(
    named_medium$evidence_source == "literature_backed_medium_catalog"
  ))
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

test_that("sample by cell type is the default inference unit", {
  expect_equal(
    eval(formals(rc_run_microcompass)$unit),
    c("sample_celltype", "metacell")
  )
  expect_equal(
    eval(formals(rc_run_regcompass)$inference_unit),
    c("sample_celltype", "metacell")
  )
  expect_identical(
    eval(formals(rc_run_regcompass_one_shot)$medium_scenario),
    "physiologic"
  )
})
