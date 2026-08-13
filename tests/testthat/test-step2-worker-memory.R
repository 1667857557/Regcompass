test_that("Step 2 worker is namespace-scoped and receives only a compact task", {
  worker <- get(
    ".rc_step2_reaction_batch_worker",
    envir = asNamespace("RegCompassR"),
    inherits = FALSE
  )
  expect_true(isNamespace(environment(worker)))
  expect_identical(names(formals(worker)), "task")
})

test_that("Step 2 compact vmax drops solver flux and unrelated attributes", {
  value <- list(
    feasible = TRUE,
    vmax = 3,
    status = "optimal",
    flux = stats::setNames(seq_len(1000), paste0("R", seq_len(1000)))
  )
  attr(value, "large_unused_attribute") <- raw(10000)
  compact <- RegCompassR:::.rc_step2_compact_vmax_value(value)
  expect_identical(
    names(compact), c("feasible", "vmax", "status", "flux")
  )
  expect_true(compact$feasible)
  expect_equal(compact$vmax, 3)
  expect_identical(compact$status, "optimal")
  expect_length(compact$flux, 0L)
  expect_null(attr(compact, "large_unused_attribute", exact = TRUE))
})

test_that("compact Step 2 payload preserves the canonical LP result", {
  skip_if_not_installed("highs")
  root <- tempfile("regcompass-step2-compact-")
  dir.create(root, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

  S <- Matrix::Matrix(
    matrix(
      c(1, -1), nrow = 1,
      dimnames = list("M", c("UP", "TARGET"))
    ),
    sparse = TRUE
  )
  lb <- c(UP = 0, TARGET = 0)
  ub <- c(UP = 10, TARGET = 10)
  model <- list(
    S = S,
    lb = lb,
    ub = ub,
    is_union_gem = TRUE,
    union_gem_cell_type = "T",
    union_gem_medium_scenario = "toy",
    union_gem_scope =
      "one_cell_type_one_medium_shared_across_conditions_and_matching_metacells",
    target_status = "ok"
  )
  model_file <- file.path(root, "union_gem.rds")
  saveRDS(model, model_file)
  checksum <- unname(tools::md5sum(model_file))
  row_id <- "cell_type=T::reaction=TARGET::direction=forward::medium=toy"
  model_cache <- stats::setNames(list(list(
    file = model_file,
    file_checksum = checksum,
    cell_type = "T",
    reaction_id = "TARGET",
    target_direction = "forward",
    medium_scenario = "toy"
  )), row_id)
  unit_celltype <- c(mc1 = "T", mc2 = "T")
  penalty_matrix <- matrix(
    c(0.25, 0.50, 0.75, 0.10),
    nrow = 2,
    dimnames = list(c("UP", "TARGET"), c("mc1", "mc2"))
  )
  vmax_result <- rc_compass_vmax_directional(
    S, lb, ub, "TARGET", direction = "forward", solver = "highs"
  )
  vmax_cache <- stats::setNames(list(vmax_result), row_id)
  payload_dir <- file.path(root, "payload")
  checkpoint_dir <- file.path(root, "checkpoint")
  dir.create(payload_dir)
  dir.create(checkpoint_dir)

  payload_file <- RegCompassR:::.rc_step2_model_payload(
    model_key = model_file,
    row_ids = row_id,
    model_cache = model_cache,
    unit_celltype = unit_celltype,
    penalties = list(penalty = penalty_matrix),
    vmax_cache = vmax_cache,
    omega = 0.95,
    solver = "highs",
    flux_threshold = 1e-8,
    payload_dir = payload_dir
  )
  payload <- readRDS(payload_file)
  expect_identical(
    names(payload),
    c(
      "schema_version", "model", "reactions", "units", "penalty",
      "entries", "vmax", "omega", "solver", "flux_threshold"
    )
  )
  expect_identical(payload$units, c("mc1", "mc2"))
  expect_length(payload$vmax[[row_id]]$flux, 0L)

  checkpoints <- RegCompassR:::.rc_step2_reaction_batch_worker(list(
    payload_file = payload_file,
    row_ids = row_id,
    checkpoint_dir = checkpoint_dir
  ))
  expect_length(checkpoints, 1L)
  observed <- readRDS(checkpoints[[1L]])
  for (one_unit in colnames(penalty_matrix)) {
    reference <- rc_compass_two_step_lp_directional(
      S, lb, ub, "TARGET", penalty_matrix[, one_unit],
      target_direction = "forward", omega = 0.95, solver = "highs"
    )
    expect_identical(observed$feasible[[one_unit]], reference$feasible)
    expect_equal(observed$vmax[[one_unit]], reference$vmax, tolerance = 1e-10)
    expect_equal(
      observed$penalty[[one_unit]], reference$penalty, tolerance = 1e-9
    )
  }
})

test_that("model-scoped reaction batches cover every target exactly once", {
  row_ids <- paste0("row_", seq_len(24L))
  model_keys <- stats::setNames(
    rep(c("model_A.rds", "model_B.rds"), each = 12L),
    row_ids
  )
  batches <- RegCompassR:::.rc_step2_model_batches(
    model_keys, workers = 6L
  )
  observed <- unlist(lapply(batches, `[[`, "row_ids"), use.names = FALSE)
  expect_length(batches, 6L)
  expect_identical(length(observed), length(unique(observed)))
  expect_setequal(observed, row_ids)
  expect_true(all(vapply(batches, function(task) {
    all(model_keys[task$row_ids] == task$model_key)
  }, logical(1))))
})
