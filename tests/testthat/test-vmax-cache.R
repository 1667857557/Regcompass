test_that("cached-vmax Step 2 matches the canonical two-step LP", {
  skip_if_not_installed("highs")
  S <- Matrix::Matrix(
    matrix(c(1, -1), nrow = 1,
           dimnames = list("M", c("UP", "TARGET"))),
    sparse = TRUE
  )
  lb <- c(UP = 0, TARGET = 0)
  ub <- c(UP = 10, TARGET = 10)
  penalties <- c(UP = 0.25, TARGET = 0.5)
  vmax <- rc_compass_vmax_directional(
    S, lb, ub, "TARGET", direction = "forward", solver = "highs"
  )
  direct <- rc_compass_two_step_lp_directional(
    S, lb, ub, "TARGET", penalties,
    target_direction = "forward", omega = 0.95, solver = "highs"
  )
  cached <- RegCompassR:::.rc_compass_step2_from_vmax_directional(
    S, lb, ub, "TARGET", penalties, vmax_result = vmax,
    target_direction = "forward", omega = 0.95, solver = "highs"
  )
  expect_identical(cached$feasible, direct$feasible)
  expect_equal(cached$vmax, direct$vmax, tolerance = 1e-10)
  expect_equal(cached$penalty, direct$penalty, tolerance = 1e-10)
  expect_equal(cached$flux, direct$flux, tolerance = 1e-10)
})

test_that("microCOMPASS engine reuses structural vmax across metacells", {
  engine <- get(
    ".rc_run_microcompass_engine",
    envir = asNamespace("RegCompassR"),
    inherits = FALSE
  )
  text <- paste(deparse(body(engine)), collapse = "\n")
  expect_match(text, ".rc_build_microcompass_vmax_cache", fixed = TRUE)
  expect_match(text, ".rc_compass_step2_from_vmax_directional", fixed = TRUE)
  expect_match(text, "vmax_reused_from_shared_cache = TRUE", fixed = TRUE)
  expect_false(grepl("rc_compass_two_step_lp_directional", text, fixed = TRUE))
})

test_that("grouped vmax results retain exact directional cache row IDs", {
  row_ids <- c(
    "reaction=R1::direction=forward::medium=toy",
    "reaction=R2::direction=reverse::medium=toy"
  )
  values <- list(
    list(feasible = TRUE, vmax = 2, status = "optimal"),
    list(feasible = TRUE, vmax = 3, status = "optimal")
  )
  grouped <- stats::setNames(
    list(stats::setNames(values, row_ids)),
    "/tmp/shared-model.rds"
  )

  answer <- RegCompassR:::.rc_flatten_microcompass_vmax_cache(
    grouped,
    rev(row_ids)
  )

  expect_identical(names(answer), rev(row_ids))
  expect_identical(answer[[row_ids[[1L]]]], values[[1L]])
  expect_error(
    RegCompassR:::.rc_flatten_microcompass_vmax_cache(
      grouped,
      c(row_ids, "reaction=missing::direction=forward::medium=toy")
    ),
    "cache is incomplete"
  )
})
