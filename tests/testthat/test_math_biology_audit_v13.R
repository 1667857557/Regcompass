
test_that("penalty LP rejects incomplete alignment", {
  S <- matrix(c(-1, 1), nrow = 1, dimnames = list("m", c("R1", "R2")))
  expect_error(
    rc_build_abs_penalty_lp(S, c(R1 = 0, R2 = 0), c(R1 = 10, R2 = 10),
                            penalties = 1, target_index = 1, target_min = 1),
    "one value per reaction"
  )
  expect_error(
    rc_build_abs_penalty_lp(S, c(R1 = 0, R2 = 0), c(R1 = 10, R2 = 10),
                            penalties = c(R1 = 1), target_index = 1, target_min = 1),
    "missing reactions"
  )
})

test_that("ranged LP rows are expanded without dropping either bound", {
  A <- Matrix::Matrix(matrix(c(1, 2), nrow = 1), sparse = TRUE)
  out <- .rc_expand_ranged_constraints(A, lhs = 1, rhs = 3)
  expect_equal(nrow(out$A), 2L)
  expect_equal(out$sense, c(">", "<"))
  expect_equal(out$bound, c(1, 3))
})

test_that("suboptimal solver status is not standardized as optimal", {
  out <- rc_standardize_lp_result("suboptimal", 0, numeric(), 0)
  expect_false(identical(out$status, "optimal"))
})

test_that("condition-specific meta-module parents preserve medium bounds", {
  skip_if_not(requireNamespace("highs", quietly = TRUE))
  S <- matrix(c(-1, -1), nrow = 1,
              dimnames = list("m", c("EX_m", "R")))
  gem <- rc_make_gem(
    S, lb = c(EX_m = -1000, R = 0), ub = c(EX_m = 0, R = 1000),
    reaction_meta = data.frame(
      reaction_id = c("EX_m", "R"), role = c("exchange", "internal"),
      role_source = "curated", stringsAsFactors = FALSE
    )
  )
  membership <- data.frame(
    sample_id = c("S1", "S2"), module_id = "M1", reaction_id = "R",
    is_core = TRUE, stringsAsFactors = FALSE
  )
  medium <- data.frame(
    medium_scenario_id = "custom", exchange_reaction_id = "EX_m",
    lb = c(-1, -10), ub = 0, available = TRUE,
    condition = c("A", "B"), stringsAsFactors = FALSE
  )
  cache <- rc_build_meta_module_gem_cache(
    gem, membership, membership, medium_scenarios = medium,
    sample_conditions = c(S1 = "A", S2 = "B"),
    solver = "highs", strict = TRUE
  )
  summary <- attr(cache, "summary")
  expect_setequal(summary$condition, c("A", "B"))
  models <- lapply(summary$file, readRDS)
  expect_setequal(vapply(models, function(x) x$lb[["EX_m"]], numeric(1)), c(-1, -10))
})

test_that("multiome penalty uses active role and GPR parameters", {
  C <- matrix(0.5, nrow = 2, ncol = 1,
              dimnames = list(c("EX", "R"), "u"))
  F <- C
  roles <- data.frame(
    reaction_id = c("EX", "R"), role = c("exchange", "internal"),
    role_source = c("curated", "curated")
  )
  gpr <- data.frame(
    reaction_id = c("EX", "R"), missing_gene_fraction = c(0, 1)
  )
  out <- rc_compute_multiome_penalty(
    C, F, gpr_diagnostics = gpr, reaction_roles = roles,
    weights = c(expr = 1, confidence = 0, missing = 0, gpr_missing = 2)
  )
  expect_equal(out$penalty["EX", "u"], 0.05)
  expect_gt(out$penalty["R", "u"], -log(0.5))
  expect_match(out$evidence_policy, "not the original COMPASS")
})
