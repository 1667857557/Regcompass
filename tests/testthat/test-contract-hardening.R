test_that("canonical GPR defaults are explicit and option-independent", {
  gpr <- list(
    complex = list(c("g1", "g2")),
    isoenzyme = list(c("g1"), c("g2"))
  )
  gene_score <- matrix(
    c(0.2, 0.9),
    ncol = 1,
    dimnames = list(c("g1", "g2"), "u1")
  )

  old <- options(RegCompassR.strict_gpr_defaults = FALSE)
  on.exit(options(old), add = TRUE)
  without_option <- rc_reaction_capacity(
    gpr,
    gene_score,
    BPPARAM = FALSE
  )

  options(RegCompassR.strict_gpr_defaults = TRUE)
  with_option <- rc_reaction_capacity(
    gpr,
    gene_score,
    BPPARAM = FALSE
  )

  expect_equal(without_option["complex", "u1"], 0.2)
  expect_equal(without_option["isoenzyme", "u1"], 0.9)
  expect_equal(with_option, without_option)
})

test_that("alternative GPR heuristics require explicit arguments", {
  gpr <- list(
    reaction = list(c("g1", "g2"), c("g3"))
  )
  gene_score <- matrix(
    c(0.2, 0.9, 0.5),
    ncol = 1,
    dimnames = list(c("g1", "g2", "g3"), "u1")
  )

  canonical <- rc_reaction_capacity(
    gpr,
    gene_score,
    BPPARAM = FALSE
  )
  heuristic <- rc_reaction_capacity(
    gpr,
    gene_score,
    promiscuity_mode = "sqrt",
    and_method = "boltzmann",
    or_method = "sum_sqrtK",
    BPPARAM = FALSE
  )

  expect_equal(canonical["reaction", "u1"], 0.5)
  expect_false(isTRUE(all.equal(heuristic, canonical)))
})

test_that("GPR-subset logCPM accepts explicit full-transcriptome library sizes", {
  counts <- Matrix::Matrix(
    matrix(
      c(10, 5, 0, 5),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("G1", "G2"), c("m1", "m2"))
    ),
    sparse = TRUE
  )

  input_scope <- .rc_metacell_logcpm(counts)
  expect_identical(
    attr(input_scope, "normalization_scope"),
    "input_matrix_library_size"
  )

  library_size <- c(m1 = 100, m2 = 200)
  observed <- .rc_metacell_logcpm(
    counts,
    library_size = library_size
  )
  expected <- log1p(
    as.matrix(counts) %*% diag(1e6 / as.numeric(library_size))
  )
  dimnames(expected) <- dimnames(counts)

  expect_equal(as.matrix(observed), expected)
  expect_identical(
    attr(observed, "normalization_scope"),
    "full_transcriptome_library_size_before_gpr_filter"
  )
})

test_that("integrated workflow validates routing inputs before delegation", {
  common <- list(
    object = NULL,
    gem = list(),
    outdir = tempfile("regcompass-routing-"),
    pfm = NULL,
    genome = NULL,
    medium_scenarios = data.frame()
  )

  expect_error(
    do.call(rc_run_regcompass, c(common, list(model_mode = "invalid"))),
    "'arg' should be one of"
  )
  expect_error(
    do.call(rc_run_regcompass, c(common, list(parallel_backend = "invalid"))),
    "'arg' should be one of"
  )
  expect_error(
    do.call(rc_run_regcompass, c(common, list(inference_unit = "invalid"))),
    "'arg' should be one of"
  )
  expect_error(
    do.call(
      rc_run_regcompass,
      c(common, list(strict_biological_defaults = NA))
    ),
    "strict_biological_defaults.*TRUE or FALSE"
  )
  expect_error(
    do.call(rc_run_regcompass, c(common, list(layer2_args = "invalid"))),
    "argument bundles must be lists: layer2_args"
  )
})
