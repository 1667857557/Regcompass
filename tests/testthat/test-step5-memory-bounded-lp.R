legacy_step2_penalty <- function(
    S, lb, ub, target_reaction, penalties, vmax_result,
    target_direction = c("forward", "reverse"), omega = 0.95,
    solver = "highs") {
  target_direction <- match.arg(target_direction)
  S <- RegCompassR:::.rc_as_dgCMatrix(S)
  reactions <- colnames(S)
  lb <- rc_align_bound(lb, reactions, default = -1000, name = "lb")
  ub <- rc_align_bound(ub, reactions, default = 1000, name = "ub")
  penalties <- RegCompassR:::.rc_compass_step2_align_penalties(
    reactions, penalties
  )
  n <- ncol(S)
  m <- nrow(S)
  target_index <- match(target_reaction, reactions)
  reaction_index <- seq_len(n)
  s_nnz_per_col <- diff(S@p)
  s_j <- rep.int(reaction_index, s_nnz_per_col)
  s_i <- S@i + 1L
  A <- Matrix::sparseMatrix(
    i = c(
      s_i,
      m + reaction_index,
      m + reaction_index,
      m + n + reaction_index,
      m + n + reaction_index,
      m + 2L * n + 1L
    ),
    j = c(
      s_j,
      reaction_index,
      n + reaction_index,
      reaction_index,
      n + reaction_index,
      target_index
    ),
    x = c(
      S@x,
      rep.int(1, n), rep.int(-1, n),
      rep.int(-1, n), rep.int(-1, n),
      if (identical(target_direction, "forward")) 1 else -1
    ),
    dims = c(m + 2L * n + 1L, 2L * n),
    giveCsparse = TRUE
  )
  answer <- rc_solve_lp(
    obj = c(rep(0, n), penalties),
    A = A,
    lhs = c(rep(0, m), rep(-Inf, 2L * n), omega * vmax_result$vmax),
    rhs = c(rep(0, m), rep(0, 2L * n), Inf),
    lb = c(lb, rep(0, n)),
    ub = c(ub, pmax(abs(lb), abs(ub))),
    solver = solver
  )
  list(
    feasible = identical(answer$status, "optimal"),
    penalty = if (identical(answer$status, "optimal")) {
      max(0, as.numeric(answer$objective))
    } else {
      NA_real_
    }
  )
}

test_that("Step 5 compiles exact directional nonnegative COMPASS LP", {
  S <- Matrix::Matrix(
    matrix(c(1, -1, 0, 1, 0, -1), nrow = 2),
    sparse = TRUE,
    dimnames = list(c("m1", "m2"), c("r1", "r2", "r3"))
  )
  lb <- c(r1 = 0, r2 = -2, r3 = 0)
  ub <- c(r1 = 10, r2 = 3, r3 = 7)
  vmax_result <- list(feasible = TRUE, vmax = 2, status = "optimal")

  observed <- RegCompassR:::.rc_compass_step2_prepare(
    S, lb, ub, "r2", vmax_result,
    target_direction = "reverse", omega = 0.95
  )
  template <- observed$template

  expect_true(observed$runnable)
  expect_identical(
    template$formulation,
    "compass_directional_nonnegative_exact_v1"
  )
  expect_identical(template$reactions, colnames(S))
  expect_identical(
    template$variable_id,
    c("r1::forward", "r2::forward", "r3::forward", "r2::reverse")
  )
  expected_A <- cbind(
    S[, c("r1", "r2", "r3"), drop = FALSE],
    -S[, "r2", drop = FALSE]
  )
  expect_equal(
    unname(as.matrix(template$A)),
    unname(as.matrix(expected_A)),
    tolerance = 0
  )
  expect_equal(template$lhs, rep(0, nrow(S)))
  expect_equal(template$rhs, rep(0, nrow(S)))
  expect_equal(template$lb, c(0, 0, 0, 1.9))
  expect_equal(template$ub, c(10, 0, 7, 2))
  expect_equal(template$n_variables, 4L)
  expect_lt(template$n_variables, 2L * ncol(S))
  expect_equal(nrow(template$A), nrow(S))
  expect_equal(template$target_variable_index, 4L)
  expect_equal(template$opposite_variable_index, 2L)
})

test_that("directional compiler preserves compulsory signed bounds", {
  S <- Matrix::Matrix(
    matrix(c(1, -1), nrow = 1),
    sparse = TRUE,
    dimnames = list("m", c("forced_forward", "forced_reverse"))
  )
  lb <- c(forced_forward = 5, forced_reverse = -9)
  ub <- c(forced_forward = 7, forced_reverse = -5)

  forward <- RegCompassR:::.rc_compass_step2_prepare(
    S, lb, ub, "forced_forward",
    list(feasible = TRUE, vmax = 7, status = "optimal"),
    target_direction = "forward", omega = 0.5
  )$template
  expect_identical(
    forward$variable_id,
    c("forced_forward::forward", "forced_reverse::reverse")
  )
  expect_equal(forward$lb, c(5, 5))
  expect_equal(forward$ub, c(7, 9))

  reverse <- RegCompassR:::.rc_compass_step2_prepare(
    S, lb, ub, "forced_reverse",
    list(feasible = TRUE, vmax = 9, status = "optimal"),
    target_direction = "reverse", omega = 0.5
  )$template
  expect_equal(reverse$lb, c(5, 5))
  expect_equal(reverse$ub, c(7, 9))
})

test_that("directional Step 2 matches the legacy signed absolute-value oracle", {
  skip_if_not_installed("highs")
  S <- Matrix::Matrix(
    matrix(c(1, -1), nrow = 1,
           dimnames = list("M", c("UP", "TARGET"))),
    sparse = TRUE
  )
  lb <- c(UP = -10, TARGET = -10)
  ub <- c(UP = 10, TARGET = 10)
  penalty_sets <- list(
    c(UP = 0.2, TARGET = 0.8),
    c(UP = 1e-6, TARGET = 1),
    c(UP = 0, TARGET = 0.5)
  )

  for (direction in c("forward", "reverse")) {
    vmax <- rc_compass_vmax_directional(
      S, lb, ub, "TARGET", direction = direction, solver = "highs"
    )
    expect_true(vmax$feasible)
    for (penalties in penalty_sets) {
      observed <- RegCompassR:::.rc_compass_step2_from_vmax_directional(
        S, lb, ub, "TARGET", penalties, vmax,
        target_direction = direction, omega = 0.95, solver = "highs"
      )
      legacy <- legacy_step2_penalty(
        S, lb, ub, "TARGET", penalties, vmax,
        target_direction = direction, omega = 0.95, solver = "highs"
      )
      expect_identical(observed$feasible, legacy$feasible)
      expect_equal(observed$penalty, legacy$penalty, tolerance = 1e-9)
    }
  }
})

test_that("Step 5 Vmax workers discard unused primal flux before collection", {
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_build_microcompass_vmax_cache_core)),
    collapse = "\n"
  )
  expect_match(implementation, "flux = numeric()", fixed = TRUE)
  expect_match(implementation, "rc_parallel_lapply", fixed = TRUE)
  expect_match(implementation, ".rc_microcompass_vmax_tasks", fixed = TRUE)
})

test_that("Step 5 memory fix preserves worker-level parallel scheduling", {
  vmax_impl <- paste(
    deparse(body(RegCompassR:::.rc_build_microcompass_vmax_cache_core)),
    collapse = "\n"
  )
  step2_impl <- paste(
    deparse(body(RegCompassR:::.rc_run_celltype_microcompass_engine_reaction_core)),
    collapse = "\n"
  )
  expect_match(vmax_impl, "workers <- .rc_microcompass_worker_count", fixed = TRUE)
  expect_match(vmax_impl, "rc_parallel_lapply", fixed = TRUE)
  expect_match(step2_impl, ".rc_step2_model_batches", fixed = TRUE)
  expect_match(step2_impl, "rc_parallel_lapply", fixed = TRUE)
})

test_that("Step-2 builder removes absolute-value and target constraint rows", {
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_compass_step2_prepare)),
    collapse = "\n"
  )
  expect_match(implementation, "direction_reaction_index", fixed = TRUE)
  expect_match(implementation, "opposite_variable", fixed = TRUE)
  expect_match(implementation, "required_flux", fixed = TRUE)
  expect_false(grepl("auxiliary_upper", implementation, fixed = TRUE))
  expect_false(grepl("2L * n_reactions + 1L", implementation, fixed = TRUE))
})
