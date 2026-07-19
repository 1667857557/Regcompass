test_that("GPR logCPM uses explicit full-transcriptome library size", {
  full <- Matrix::Matrix(
    matrix(
      c(10, 10, 0, 0, 90, 990),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("G1", "G2", "OTHER"), c("m1", "m2"))
    ),
    sparse = TRUE
  )
  observed <- .rc_metacell_logcpm(
    full[c("G1", "G2"), , drop = FALSE],
    library_size = Matrix::colSums(full)
  )
  expected <- log1p(c(10 / 100, 10 / 1000) * 1e6)
  expect_equal(as.numeric(observed["G1", ]), expected)
  expect_identical(
    attr(observed, "normalization_scope"),
    "full_transcriptome_library_size_before_gpr_filter"
  )
})

test_that("missing and no-GPR expression are not cheaper than observed zero", {
  reactions <- c(
    "observed_zero", "assay_missing", "no_gpr",
    "high_expression", "exchange", "demand"
  )
  expression <- matrix(
    c(0, NA, NA, 0.8, NA, NA),
    ncol = 1,
    dimnames = list(reactions, "u1")
  )
  roles <- data.frame(
    reaction_id = reactions,
    role = c(
      "internal", "internal", "internal",
      "internal", "exchange", "demand"
    ),
    role_source = c(
      "metadata", "metadata", "unknown",
      "metadata", "metadata", "id_pattern"
    ),
    stringsAsFactors = FALSE
  )
  answer <- rc_compute_multiome_penalty(expression, reaction_roles = roles)

  expect_equal(answer$penalty["observed_zero", "u1"], 1)
  expect_equal(answer$penalty["assay_missing", "u1"], 1)
  expect_equal(answer$penalty["no_gpr", "u1"], 1)
  expect_equal(
    answer$penalty["high_expression", "u1"],
    1 / (1 + log2(1 + 0.8))
  )
  expect_equal(answer$penalty["exchange", "u1"], 1)
  expect_equal(answer$penalty["demand", "u1"], 20)
})

test_that("species-matched physiological media remain defaults", {
  expect_identical(eval(formals(rc_make_medium_scenarios)$scenario),
                   "physiologic")
  expect_identical(eval(formals(rc_run_regcompass_one_shot)$medium_scenario),
                   "physiologic")
})

test_that("main workflow routes Pando by condition and cell type", {
  body_text <- paste(deparse(body(rc_run_regcompass)), collapse = "\n")
  expect_match(body_text, ".rc_run_condition_pando_modules", fixed = TRUE)
  expect_match(body_text, "condition_col", fixed = TRUE)
  expect_match(body_text, "celltype_col", fixed = TRUE)
})
