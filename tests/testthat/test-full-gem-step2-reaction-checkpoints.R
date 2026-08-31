test_that("full-GEM Step 2 uses reaction batches and file-backed checkpoints", {
  core <- paste(
    deparse(body(RegCompassR:::.rc_run_shared_full_gem_engine_core)),
    collapse = "\n"
  )
  worker <- RegCompassR:::.rc_full_gem_step2_reaction_batch_worker

  expect_true(isNamespace(environment(worker)))
  expect_identical(names(formals(worker)), "task")
  expect_match(core, ".rc_step2_model_batches", fixed = TRUE)
  expect_match(
    core, ".rc_full_gem_step2_reaction_batch_worker", fixed = TRUE
  )
  expect_match(core, ".rc_atomic_save_rds", fixed = TRUE)
  expect_match(core, "unlink(checkpoint_files[[i]]", fixed = TRUE)
  expect_false(grepl(
    "expand.grid(\\n    model_key = unique_model_keys,\\n    unit_id = units",
    core, fixed = TRUE
  ))
})

test_that("full-GEM reaction worker preserves canonical Step 2 LP values", {
  skip_if_not_installed("highs")

  root <- tempfile("regcompass-full-gem-step2-")
  dir.create(root, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

  S <- Matrix::Matrix(
    matrix(
      c(1, -1, -1),
      nrow = 1,
      dimnames = list("M", c("UP", "T1", "T2"))
    ),
    sparse = TRUE
  )
  lb <- c(UP = 0, T1 = 0, T2 = 0)
  ub <- c(UP = 10, T1 = 10, T2 = 10)
  model <- list(
    S = S,
    lb = lb,
    ub = ub,
    target_status = "not_prechecked",
    closure_diagnostics = data.frame()
  )
  model_file <- file.path(root, "full_gem.rds")
  saveRDS(model, model_file)

  row_ids <- c(
    "reaction=T1::direction=forward::medium=toy::condition=all",
    "reaction=T2::direction=forward::medium=toy::condition=all"
  )
  model_cache <- stats::setNames(lapply(c("T1", "T2"), function(target) {
    list(
      reaction_id = target,
      target_direction = "forward",
      medium_scenario = "toy",
      condition = "all",
      file = model_file,
      build_strategy = "compass_medium_constrained_full_gem"
    )
  }), row_ids)

  units <- c("mc1", "mc2", "mc3")
  penalty_matrix <- matrix(
    c(
      0.20, 0.50, 0.90,
      0.25, 0.80, 0.40,
      0.70, 0.10, 0.60
    ),
    nrow = 3,
    dimnames = list(c("UP", "T1", "T2"), units)
  )
  vmax_cache <- stats::setNames(lapply(c("T1", "T2"), function(target) {
    rc_compass_vmax_directional(
      S = S,
      lb = lb,
      ub = ub,
      target_reaction = target,
      direction = "forward",
      solver = "highs"
    )
  }), row_ids)

  payload_dir <- file.path(root, "payload")
  checkpoint_dir <- file.path(root, "checkpoint")
  dir.create(payload_dir)
  dir.create(checkpoint_dir)

  payload_file <- RegCompassR:::.rc_full_gem_step2_model_payload(
    model_key = model_file,
    row_ids = row_ids,
    model_cache = model_cache,
    penalties = list(penalty = penalty_matrix),
    vmax_cache = vmax_cache,
    omega = 0.95,
    solver = "highs",
    flux_threshold = 1e-8,
    payload_dir = payload_dir
  )
  payload <- readRDS(payload_file)
  expect_identical(
    payload$schema_version,
    "regcompass_full_gem_step2_compact_payload_v1"
  )
  expect_identical(payload$units, units)
  expect_identical(payload$penalty, penalty_matrix)
  expect_true(all(vapply(payload$vmax, function(x) {
    length(x$flux) == 0L
  }, logical(1))))

  checkpoints <- RegCompassR:::.rc_full_gem_step2_reaction_batch_worker(
    list(
      payload_file = payload_file,
      row_ids = row_ids,
      checkpoint_dir = checkpoint_dir
    )
  )
  expect_length(checkpoints, length(row_ids))
  expect_true(all(file.exists(checkpoints)))

  observed_rows <- character()
  for (checkpoint in checkpoints) {
    observed <- readRDS(checkpoint)
    row_id <- observed$row_id
    observed_rows <- c(observed_rows, row_id)
    target <- model_cache[[row_id]]$reaction_id

    expect_identical(observed$units, units)
    expect_equal(observed$engine_metrics$n_solves, length(units))
    expect_equal(nrow(observed$diagnostics), length(units))

    for (one_unit in units) {
      reference <- RegCompassR:::.rc_compass_step2_from_vmax_directional(
        S = S,
        lb = lb,
        ub = ub,
        target_reaction = target,
        penalties = penalty_matrix[, one_unit],
        vmax_result = vmax_cache[[row_id]],
        target_direction = "forward",
        omega = 0.95,
        solver = "highs",
        flux_threshold = 1e-8
      )
      expect_identical(
        observed$feasible[[one_unit]], reference$feasible
      )
      expect_equal(
        observed$vmax[[one_unit]], reference$vmax,
        tolerance = 1e-10
      )
      expect_equal(
        observed$penalty[[one_unit]], reference$penalty,
        tolerance = 1e-9
      )
    }
  }
  expect_setequal(observed_rows, row_ids)
})

test_that("full-GEM Step 2 batching preserves worker-width reaction parallelism", {
  row_ids <- paste0("row_", seq_len(200L))
  model_keys <- stats::setNames(rep("full_gem.rds", length(row_ids)), row_ids)
  batches <- RegCompassR:::.rc_step2_model_batches(
    model_keys, workers = 80L
  )
  observed <- unlist(lapply(batches, `[[`, "row_ids"), use.names = FALSE)

  expect_length(batches, 80L)
  expect_identical(length(observed), length(unique(observed)))
  expect_setequal(observed, row_ids)
})
