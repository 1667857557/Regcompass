
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

test_that("suboptimal solver status is not classified as optimal", {
  expect_false(identical(.rc_lp_status("suboptimal"), "optimal"))
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
  expect_equal(out$penalty["EX", "u"], 1)
  expect_gt(out$penalty["R", "u"], 1 - 0.5)
  expect_match(out$evidence_description, "not the original COMPASS")
})

test_that("sample aggregation preserves continuous covariate classes", {
  score <- matrix(1:8, nrow = 1,
                  dimnames = list("reaction=R::direction=forward::medium=base", paste0("u", 1:8)))
  meta <- data.frame(
    unit_id = paste0("u", 1:8),
    sample_id = rep(c("S1", "S2", "S3", "S4"), each = 2),
    condition = rep(c("A", "A", "B", "B"), each = 2),
    age = rep(c(30, 40, 50, 60), each = 2),
    stringsAsFactors = FALSE
  )
  out <- .rc_aggregate_microcompass_samples(
    score, meta, "sample_id", "condition", covariates = "age"
  )
  expect_true(is.numeric(out$meta$age))
  expect_equal(out$meta$age, c(30, 40, 50, 60))
})

test_that("blocked full-GEM directions are retained only as diagnostics", {
  S <- matrix(0, nrow = 1, ncol = 2,
              dimnames = list("m", c("blocked", "open")))
  gem <- rc_make_gem(S, lb = c(blocked = 0, open = 0),
                     ub = c(blocked = 0, open = 10))
  directions <- rc_prepare_directional_targets(
    gem, c("blocked", "open"), target_direction = "both"
  )
  allowed <- directions[directions$target_direction %in% c("forward", "reverse"), , drop = FALSE]
  expect_equal(directions$target_direction[directions$reaction_id == "blocked"], "none")
  expect_false("blocked" %in% allowed$reaction_id)
  expect_equal(allowed$target_direction, "forward")
})
