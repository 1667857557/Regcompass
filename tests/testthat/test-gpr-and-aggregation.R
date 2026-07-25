test_that("GPR-AND supports the three COMPASS functions", {
  scores <- c(0.2, 0.5, 0.9)
  expect_equal(rc_and_capacity(scores, method = "min"), 0.2)
  expect_equal(rc_and_capacity(scores, method = "median"), 0.5)
  expect_equal(rc_and_capacity(scores, method = "mean"), mean(scores))
})

test_that("RegCompass defaults GPR-AND aggregation to min", {
  parsed <- list(complex = list(c("g1", "g2", "g3")))
  gene_score <- matrix(
    c(0.2, 0.5, 0.9),
    ncol = 1,
    dimnames = list(c("g1", "g2", "g3"), "u1")
  )
  observed <- rc_reaction_capacity(
    parsed,
    gene_score,
    promiscuity_mode = "none",
    or_method = "sum",
    BPPARAM = FALSE
  )
  expect_equal(observed["complex", "u1"], 0.2)
  expect_identical(eval(formals(rc_and_capacity)$method), c("min", "median", "mean"))
  expect_identical(
    eval(formals(rc_reaction_capacity)$and_method),
    c("min", "median", "mean")
  )
  expect_identical(
    eval(formals(rc_regcompass_step_layer1)$gpr_and_method),
    c("min", "median", "mean")
  )
})

test_that("median and mean are propagated through reaction GPR evaluation", {
  parsed <- list(c("g1", "g2", "g3"), "g4")
  gene_score <- c(g1 = 0.2, g2 = 0.5, g3 = 0.9, g4 = 0.4)

  median_capacity <- rc_reaction_capacity_one(
    parsed,
    gene_score,
    and_method = "median",
    or_method = "sum"
  )
  mean_capacity <- rc_reaction_capacity_one(
    parsed,
    gene_score,
    and_method = "mean",
    or_method = "sum"
  )

  expect_equal(median_capacity, stats::median(c(0.2, 0.5, 0.9)) + 0.4)
  expect_equal(mean_capacity, mean(c(0.2, 0.5, 0.9)) + 0.4)
})

test_that("retired Boltzmann and tau APIs are absent", {
  expect_false(exists(
    "rc_boltzmann_minavg",
    envir = asNamespace("RegCompassR"),
    inherits = FALSE
  ))
  expect_false("tau" %in% names(formals(rc_reaction_capacity)))
  expect_false("tau" %in% names(formals(rc_reaction_capacity_one)))
  expect_false("tau" %in% names(formals(rc_and_capacity)))
  expect_false("tau" %in% names(formals(rc_regcompass_step_layer1)))
  expect_error(
    rc_and_capacity(c(0.2, 0.8), method = "boltzmann"),
    "should be one of"
  )
})

test_that("reaction evidence recomputation uses the stored AND method", {
  layer1 <- list(
    parsed_gpr = list(R1 = list(c("g1", "g2", "g3"))),
    gene_support_rna = matrix(
      c(0.2, 0.5, 0.9),
      ncol = 1,
      dimnames = list(c("g1", "g2", "g3"), "u1")
    ),
    gene_support_multiome = matrix(
      c(0.3, 0.6, 1),
      ncol = 1,
      dimnames = list(c("g1", "g2", "g3"), "u1")
    ),
    capacity_params = list(
      promiscuity_mode = "none",
      and_method = "median",
      or_method = "sum"
    )
  )
  catalog <- data.frame(reaction_id = "R1", stringsAsFactors = FALSE)
  capacities <- .rc_ra_reaction_capacity_pair(catalog, layer1)
  expect_equal(capacities$rna["R1", "u1"], 0.5)
  expect_equal(capacities$multiome["R1", "u1"], 0.6)
})
