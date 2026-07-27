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

test_that("COMPASS GPR-AND alternatives require explicit arguments", {
  gpr <- list(
    reaction = list(c("g1", "g2", "g3"))
  )
  gene_score <- matrix(
    c(0.2, 0.5, 0.9),
    ncol = 1,
    dimnames = list(c("g1", "g2", "g3"), "u1")
  )

  minimum <- rc_reaction_capacity(
    gpr,
    gene_score,
    promiscuity_mode = "none",
    or_method = "sum",
    BPPARAM = FALSE
  )
  median <- rc_reaction_capacity(
    gpr,
    gene_score,
    promiscuity_mode = "none",
    and_method = "median",
    or_method = "sum",
    BPPARAM = FALSE
  )
  average <- rc_reaction_capacity(
    gpr,
    gene_score,
    promiscuity_mode = "none",
    and_method = "mean",
    or_method = "sum",
    BPPARAM = FALSE
  )

  expect_equal(minimum["reaction", "u1"], 0.2)
  expect_equal(median["reaction", "u1"], 0.5)
  expect_equal(average["reaction", "u1"], mean(c(0.2, 0.5, 0.9)))
  expect_error(
    rc_reaction_capacity(
      gpr,
      gene_score,
      and_method = "boltzmann",
      BPPARAM = FALSE
    ),
    "should be one of"
  )
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

  triangular <- Matrix::triu(Matrix::Matrix(
    matrix(c(1, 1, 0, 1), nrow = 2),
    sparse = TRUE
  ))
  expect_silent(triangular_output <- .rc_metacell_logcpm(triangular))
  expect_equal(
    as.matrix(triangular_output),
    log1p(as.matrix(triangular) %*%
            diag(1e6 / Matrix::colSums(triangular)))
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
    do.call(rc_run_regcompass, c(common, list(layer2_args = "invalid"))),
    "argument bundles must be lists: layer2_args"
  )
  expect_error(
    do.call(
      rc_run_regcompass,
      c(common, list(layer1_args = list(tau = 0.2)))
    ),
    "Unknown `layer1_args` fields: tau"
  )
  expect_error(.rc_stage_worker_config(0L, "upstream_workers"), "at least 1")
  expect_error(.rc_stage_worker_config(0L, "layer2_workers"), "at least 1")

  formals_names <- names(formals(rc_run_regcompass))
  expect_true("model_mode" %in% formals_names)
  expect_true("upstream_workers" %in% formals_names)
  expect_true("layer2_workers" %in% formals_names)
  expect_false("parallel_backend" %in% formals_names)
  expect_false("inference_unit" %in% formals_names)
  expect_false("strict_biological_defaults" %in% formals_names)
  expect_identical(formals(rc_run_regcompass)$upstream_workers, 6L)
  expect_identical(formals(rc_run_regcompass)$layer2_workers, 30L)
  expect_identical(
    eval(formals(rc_run_regcompass)$model_mode),
    c("meta_module_gem", "full_gem")
  )
})
