rc_test_toy_gem <- function() {
  S <- matrix(
    c(
      1, -1, 0,
      0,  1, -1
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("m1", "m2"),
      c("EX_m1", "Rtarget", "EX_m2")
    )
  )
  rc_make_gem(
    S,
    lb = c(EX_m1 = 0, Rtarget = 0, EX_m2 = 0),
    ub = c(EX_m1 = 1000, Rtarget = 1000, EX_m2 = 1000),
    reaction_meta = data.frame(
      reaction_id = colnames(S),
      role = c("exchange", "internal", "exchange"),
      role_source = "curated",
      stringsAsFactors = FALSE
    )
  )
}

test_that("reverse vmax uses the minimization objective convention", {
  S <- matrix(
    0,
    nrow = 1,
    ncol = 1,
    dimnames = list("dummy", "Rrev")
  )
  answer <- rc_compass_vmax_directional(
    S,
    lb = c(Rrev = -3),
    ub = c(Rrev = 2),
    target_reaction = "Rrev",
    direction = "reverse"
  )
  if (identical(answer$status, "error")) {
    skip("No LP solver available")
  }
  expect_true(answer$feasible)
  expect_equal(answer$vmax, 3, tolerance = 1e-6)
})

test_that("curated support penalty overrides evidence penalty", {
  capacity <- matrix(
    0.01,
    nrow = 3,
    ncol = 1,
    dimnames = list(
      c("EX_curated", "TR_curated", "EX_inferred"),
      "u1"
    )
  )
  confidence <- capacity
  roles <- data.frame(
    reaction_id = rownames(capacity),
    role = c("exchange", "transport", "exchange"),
    role_source = c("curated", "curated", "id_pattern"),
    stringsAsFactors = FALSE
  )
  output <- rc_compute_multiome_penalty(
    capacity,
    confidence,
    reaction_roles = roles
  )
  expect_equal(output$evidence_policy, "penalty_only")
  expect_equal(output$penalty["EX_curated", "u1"], 0.05)
  expect_gt(output$penalty["EX_inferred", "u1"], 0.05)
})

test_that("single-stoichiometry reactions are boundary-like", {
  S <- matrix(
    c(1, 0, 1, -1),
    nrow = 2,
    dimnames = list(
      c("m1", "m2"),
      c("R_boundary", "R_internal")
    )
  )
  gem <- rc_annotate_reaction_roles(
    rc_make_gem(
      S,
      lb = c(0, 0),
      ub = c(1000, 1000)
    ),
    infer_from_id = FALSE,
    infer_from_compartment = FALSE
  )
  roles <- stats::setNames(
    gem$reaction_roles$role,
    gem$reaction_roles$reaction_id
  )
  expect_equal(roles[["R_boundary"]], "boundary_like")
})

test_that("row IDs require and parse v1.3 labeled format", {
  expect_error(
    rc_parse_microcompass_row_id("R1::forward::blood_like"),
    "v1.3 labeled format"
  )
  expect_error(
    rc_parse_microcompass_row_id("reaction=R1::direction=both::medium=base"),
    "direction"
  )
  expect_error(
    rc_parse_microcompass_row_id("reaction=::direction=forward::medium=base"),
    "non-empty"
  )
  expect_error(
    rc_parse_microcompass_row_id(paste0(
      "reaction=R1::reaction=R2::direction=forward::medium=base"
    )),
    "exactly one"
  )
  parsed <- rc_parse_microcompass_row_id(c(
    paste0(
      "reaction=R1::direction=forward",
      "::medium=blood_like::condition=ctrl"
    ),
    paste0(
      "sample=S1::module=S1%3A%3AM1::reaction=R3",
      "::direction=reverse::medium=base::condition=treat"
    )
  ))
  expect_equal(parsed$reaction_id, c("R1", "R3"))
  expect_equal(
    parsed$target_direction,
    c("forward", "reverse")
  )
  expect_equal(
    parsed$medium_scenario,
    c("blood_like", "base")
  )
  expect_equal(parsed$condition, c("ctrl", "treat"))
  expect_equal(parsed$sample_id[[2L]], "S1")
  expect_equal(parsed$module_id[[2L]], "S1::M1")
})

test_that("full-GEM cache is structural and evidence-independent", {
  gem <- rc_test_toy_gem()
  directions <- data.frame(
    reaction_id = "Rtarget",
    target_direction = "forward",
    stringsAsFactors = FALSE
  )
  first <- rc_build_full_gem_cache(gem, directions, NULL)
  second <- rc_build_full_gem_cache(gem, directions, NULL)
  expect_identical(
    lapply(first, function(entry) colnames(readRDS(entry$file)$S)),
    lapply(second, function(entry) colnames(readRDS(entry$file)$S))
  )
})

test_that("penalty is unit-specific", {
  capacity <- matrix(
    c(0.9, 0.1),
    nrow = 1,
    dimnames = list("Rtarget", c("u1", "u2"))
  )
  confidence <- matrix(
    1,
    nrow = 1,
    ncol = 2,
    dimnames = dimnames(capacity)
  )
  penalty <- rc_compute_multiome_penalty(
    capacity,
    confidence
  )
  expect_false(identical(
    penalty$penalty[, 1L],
    penalty$penalty[, 2L]
  ))
})

test_that("parallel and serial full-GEM scoring agree", {
  gem <- rc_test_toy_gem()
  layer1 <- list(
    C_rel = matrix(
      c(0.9, 0.2, 0.8, 0.8, 0.7, 0.3),
      nrow = 3,
      dimnames = list(
        c("EX_m1", "Rtarget", "EX_m2"),
        c("p1", "p2")
      )
    ),
    reaction_confidence = matrix(
      1,
      nrow = 3,
      ncol = 2,
      dimnames = list(
        c("EX_m1", "Rtarget", "EX_m2"),
        c("p1", "p2")
      )
    ),
    gpr_diagnostics = NULL,
    unit_meta = data.frame(
      pool_id = c("p1", "p2"),
      sample_id = c("s1", "s2"),
      condition = c("a", "b"),
      cell_type = "T",
      stringsAsFactors = FALSE
    )
  )
  serial <- rc_run_microcompass(
    layer1,
    gem,
    "Rtarget",
    mode = "full_gem",
    unit = "metacell",
    target_direction = "forward",
    parallel = FALSE
  )
  parallel <- rc_run_microcompass(
    layer1,
    gem,
    "Rtarget",
    mode = "full_gem",
    unit = "metacell",
    target_direction = "forward",
    parallel = TRUE,
    BPPARAM = FALSE
  )
  expect_equal(serial$score, parallel$score, tolerance = 1e-8)
  expect_equal(serial$penalty, parallel$penalty, tolerance = 1e-8)
  expect_equal(serial$feasible, parallel$feasible)
})
