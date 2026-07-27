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
  out <- .rc_integrate_regulatory_support(C, R, alpha = 1)
  expect_equal(out["zero", "u1"], 0)
  expect_gt(out["activated", "u1"], C["activated", "u1"])
  expect_lt(out["repressed", "u1"], C["repressed", "u1"])
  expect_true(all(out >= 0 & out <= 1))
})

test_that("COMPASS GPR-AND functions are exact", {
  x <- c(0.2, 0.5, 0.9)
  expect_equal(rc_and_capacity(x), min(x))
  expect_equal(rc_and_capacity(x, "median"), stats::median(x))
  expect_equal(rc_and_capacity(x, "mean"), mean(x))
  expect_error(rc_and_capacity(x, "boltzmann"), "should be one of")
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
  expect_equal(
    answer$components$effective_reaction_expression["missing", "u1"], 0
  )
  expect_gt(penalty[["zero"]], penalty[["low"]])
  expect_gt(penalty[["low"]], penalty[["high"]])
  expect_true(all(is.finite(penalty) & penalty > 0))
  expect_identical(answer$evidence_policy, "penalty_only")
  expect_identical(
    answer$penalty_version,
    "gene_integrated_multiome_penalty_v1"
  )
})

test_that("Stage 1 installs species-specific Pando motifs and regions", {
  text <- paste(
    deparse(body(.rc_run_condition_single_cell_grns_legacy)),
    collapse = "\n"
  )
  motif_helper <- paste(
    deparse(body(.rc_default_pando_motifs)), collapse = "\n"
  )
  region_helper <- paste(
    deparse(body(.rc_default_pando_regions)), collapse = "\n"
  )
  expect_match(text, ".rc_default_pando_motifs", fixed = TRUE)
  expect_match(text, ".rc_default_pando_regions(species)", fixed = TRUE)
  expect_match(motif_helper, 'list = "motifs"', fixed = TRUE)
  expect_match(region_helper, "phastConsElements20Mammals.UCSC.hg38", fixed = TRUE)
  expect_match(region_helper, "SCREEN.ccRE.UCSC.hg38", fixed = TRUE)
  expect_match(region_helper, 'identical(species, "mouse")', fixed = TRUE)
  expect_match(region_helper, "return(phast_cons)", fixed = TRUE)
  expect_match(region_helper, "BiocGenerics::union", fixed = TRUE)
})

test_that("canonical scoring is metacell only", {
  expect_false("inference_unit" %in% names(formals(rc_run_regcompass)))
  expect_identical(eval(formals(rc_run_microcompass)$unit), "metacell")
  expect_false("sample_col" %in% names(formals(rc_run_microcompass)))
  expect_error(
    rc_run_microcompass(NULL, NULL, unit = "sample_celltype"),
    "historical sample-by-cell-type aggregation mode"
  )
  expect_identical(
    eval(formals(rc_run_regcompass_one_shot)$medium_scenario),
    "physiologic"
  )
})
