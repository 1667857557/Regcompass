test_that("native CORDA2 flux scanning matches ordered R semantics", {
  flux <- c(0, 2, 3, 1e-9, 4, NA_real_)
  class_code <- c(1L, 2L, 3L, 4L, 2L, 3L)
  track_code <- c(2L, 3L)
  threshold <- 1e-8

  native <- RegCompassR:::.rc_corda2_scan_flux_cpp(
    flux, class_code, track_code, threshold
  )
  expected_active <- which(is.finite(flux) & flux > threshold)
  expected_used <- expected_active[
    class_code[expected_active] %in% track_code
  ]

  expect_identical(as.integer(native$active), as.integer(expected_active))
  expect_identical(as.integer(native$used), as.integer(expected_used))
})

test_that("vmax batching fills workers within a shared structural model", {
  row_ids <- paste0("row_", seq_len(12L))
  model_keys <- stats::setNames(
    rep("/tmp/one-shared-model.rds", length(row_ids)), row_ids
  )

  tasks <- RegCompassR:::.rc_microcompass_vmax_tasks(
    model_keys, workers = 4L
  )
  observed <- unlist(lapply(tasks, `[[`, "row_ids"), use.names = FALSE)

  expect_equal(length(tasks), 4L)
  expect_setequal(observed, row_ids)
  expect_equal(length(observed), length(unique(observed)))
  expect_true(all(vapply(
    tasks,
    function(task) identical(task$model_key, "/tmp/one-shared-model.rds"),
    logical(1)
  )))
})

test_that("persistent Step 2 reuses one native model without score drift", {
  skip_if_not_installed("highs")
  S <- Matrix::Matrix(
    matrix(c(1, -1), nrow = 1,
           dimnames = list("M", c("UP", "TARGET"))),
    sparse = TRUE
  )
  lb <- c(UP = 0, TARGET = 0)
  ub <- c(UP = 10, TARGET = 10)
  vmax <- rc_compass_vmax_directional(
    S, lb, ub, "TARGET", direction = "forward", solver = "highs"
  )
  prepared <- RegCompassR:::.rc_compass_step2_prepare(
    S, lb, ub, "TARGET", vmax,
    target_direction = "forward", omega = 0.95
  )
  engine <- RegCompassR:::.rc_compass_step2_new_engine(
    prepared$template, "highs"
  )
  on.exit(
    RegCompassR:::.rc_compass_step2_release_engine(engine),
    add = TRUE
  )
  if (!identical(engine$type, "highs_persistent_cpp")) {
    skip("Installed highs does not expose the persistent native API")
  }

  penalty_sets <- list(
    c(UP = 0.25, TARGET = 0.50),
    c(UP = 0.75, TARGET = 0.10),
    c(UP = 0.05, TARGET = 1.25)
  )
  for (penalties in penalty_sets) {
    solved <- RegCompassR:::.rc_compass_step2_engine_solve(
      engine, penalties
    )
    engine <- solved$engine
    persistent <- RegCompassR:::.rc_compass_step2_result(
      prepared$template, solved$answer
    )
    reference <- rc_compass_two_step_lp_directional(
      S, lb, ub, "TARGET", penalties,
      target_direction = "forward", omega = 0.95, solver = "highs"
    )

    expect_identical(persistent$feasible, reference$feasible)
    expect_equal(persistent$vmax, reference$vmax, tolerance = 1e-10)
    expect_equal(persistent$penalty, reference$penalty, tolerance = 1e-9)
  }

  metrics <- RegCompassR:::.rc_compass_step2_engine_metrics(engine)
  expect_identical(metrics$engine, "highs_persistent_cpp")
  expect_equal(metrics$n_solves, length(penalty_sets))
  expect_gt(metrics$n_objective_updates, 0L)
  expect_equal(metrics$n_fallback, 0L)
})
