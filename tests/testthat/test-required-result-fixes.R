test_that("GPR logCPM uses full-transcriptome library size", {
  full <- Matrix::Matrix(
    matrix(
      c(
        10, 10,
        0, 0,
        90, 990
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("G1", "G2", "OTHER"), c("m1", "m2"))
    ),
    sparse = TRUE
  )
  subset_counts <- full[c("G1", "G2"), , drop = FALSE]
  key <- .rc_library_cache_key(colnames(full))
  assign(
    key,
    Matrix::colSums(full),
    envir = .rc_full_library_size_cache
  )

  observed <- .rc_metacell_logcpm(subset_counts)
  expected <- log1p(c(10 / 100, 10 / 1000) * 1e6)

  expect_equal(as.numeric(observed["G1", ]), expected)
  expect_identical(
    attr(observed, "normalization_scope"),
    "full_transcriptome_library_size_before_gpr_filter"
  )
  expect_false(exists(
    key,
    envir = .rc_full_library_size_cache,
    inherits = FALSE
  ))
})

test_that("missing and no-GPR evidence are not cheaper than observed zero", {
  reactions <- c(
    "observed_zero", "assay_missing", "no_gpr",
    "high_support", "exchange", "demand"
  )
  capacity <- matrix(
    c(0, NA, NA, 0.8, NA, NA),
    ncol = 1,
    dimnames = list(reactions, "u1")
  )
  regulation <- matrix(
    0.5,
    nrow = length(reactions),
    ncol = 1,
    dimnames = dimnames(capacity)
  )
  diagnostics <- data.frame(
    reaction_id = c("observed_zero", "assay_missing", "high_support"),
    missing_gene_fraction = c(0, 1, 0),
    stringsAsFactors = FALSE
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
    role_confidence = c(
      "medium", "medium", "low",
      "medium", "medium", "low"
    ),
    stringsAsFactors = FALSE
  )

  answer <- rc_compute_multiome_penalty(
    capacity,
    regulation,
    gpr_diagnostics = diagnostics,
    reaction_roles = roles
  )

  expect_equal(answer$penalty["observed_zero", "u1"], 1)
  expect_equal(answer$penalty["assay_missing", "u1"], 1)
  expect_equal(answer$penalty["no_gpr", "u1"], 1)
  expect_equal(answer$penalty["high_support", "u1"], 0.2)
  expect_equal(answer$penalty["exchange", "u1"], 1)
  expect_equal(answer$penalty["demand", "u1"], 20)
  expect_gte(
    answer$penalty["assay_missing", "u1"],
    answer$penalty["observed_zero", "u1"]
  )
  expect_gte(
    answer$penalty["no_gpr", "u1"],
    answer$penalty["observed_zero", "u1"]
  )
})

test_that("COMPASS-style medium preserves model direction and caps exchanges", {
  S <- Matrix::Matrix(
    matrix(c(-1, -1, -1, 1), nrow = 1),
    sparse = TRUE
  )
  rownames(S) <- "m_e"
  colnames(S) <- c("EX_both", "EX_secrete", "EX_uptake", "R1")
  gem <- list(
    S = S,
    lb = stats::setNames(c(-1000, 0, -1000, 0), colnames(S)),
    ub = stats::setNames(c(1000, 1000, 0, 1000), colnames(S)),
    reaction_meta = data.frame(
      reaction_id = colnames(S),
      role = c("exchange", "exchange", "exchange", "internal"),
      stringsAsFactors = FALSE
    )
  )

  medium <- rc_make_medium_scenarios(
    gem,
    scenario = "compass_model_bounds"
  )
  medium <- medium[
    match(c("EX_both", "EX_secrete", "EX_uptake"),
          medium$exchange_reaction_id),
    ,
    drop = FALSE
  ]

  expect_equal(medium$lb, c(-1, 0, -1))
  expect_equal(medium$ub, c(1, 1, 0))
  expect_true(all(medium$condition == "all"))
  expect_true(all(medium$exchange_limit == 1))
  expect_true(all(
    medium$assumption_level == "shared_model_defined_environment"
  ))
})

test_that("COMPASS-style model bounds are the workflow defaults", {
  expect_identical(
    eval(formals(rc_make_medium_scenarios)$scenario),
    "compass_model_bounds"
  )
  expect_identical(
    eval(formals(rc_run_regcompass_one_shot)$medium_scenario),
    "compass_model_bounds"
  )
})

test_that("peak-gene inference remains strict-stratum specific", {
  body_text <- paste(deparse(body(rc_run_pando_meta_modules)), collapse = "\n")
  expect_match(body_text, "\\.rc_pando_group_id")
  expect_match(body_text, "run_one_group")
  expect_match(body_text, "subset\\(metacell_object, cells = cells\\)")
})
