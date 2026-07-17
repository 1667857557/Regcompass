test_that("workflow OR aggregation matches the recorded parameter", {
  expect_equal(
    .rc_resolve_workflow_or_method(
      list(),
      strict_biological_defaults = TRUE
    ),
    "max"
  )
  expect_equal(
    .rc_workflow_or_method_source(
      list(),
      strict_biological_defaults = TRUE
    ),
    "strict_biological_default"
  )
  expect_equal(
    .rc_resolve_workflow_or_method(list(
      promiscuity_mode = "none",
      and_method = "min",
      or_method = "sum_sqrtK"
    )),
    "sum_sqrtK"
  )

  artifact_file <- tempfile(fileext = ".rds")
  saveRDS(
    list(capacity_params = list(or_method = "sum_sqrtK")),
    artifact_file
  )
  used <- .rc_finalize_stratum_capacity_params(
    artifact_file,
    list(promiscuity_mode = "none", and_method = "min")
  )
  artifact <- readRDS(artifact_file)

  expect_equal(used, "max")
  expect_equal(artifact$capacity_params$or_method, "max")
  expect_equal(
    artifact$capacity_params$or_method_source,
    "strict_biological_default"
  )
})


test_that("Q95 shrinkage is diagnostic and its arguments are effective", {
  C_raw <- rbind(
    r1 = c(0.10, 0.20, 0.80, 0.90),
    r2 = c(0.20, 0.30, 0.60, 0.70)
  )
  colnames(C_raw) <- paste0("p", 1:4)
  unit_meta <- data.frame(
    pool_id = colnames(C_raw),
    cell_type = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  unshrunk <- rc_q95_calibrate(
    C_raw,
    bootstrap = FALSE,
    n0 = 0,
    unit_meta = unit_meta,
    stratum_col = "cell_type"
  )
  shrunk <- rc_q95_calibrate(
    C_raw,
    bootstrap = FALSE,
    n0 = 1000,
    unit_meta = unit_meta,
    stratum_col = "cell_type"
  )

  expect_equal(unshrunk$C_rel, C_raw)
  expect_equal(shrunk$C_rel, C_raw)
  expect_false(isTRUE(all.equal(
    unshrunk$C_within_reaction_relative,
    shrunk$C_within_reaction_relative
  )))
  expect_true(all(
    shrunk$Q$calibration_role == "diagnostic_only_not_lp_capacity"
  ))
  expect_true(all(shrunk$Q$rho_n < unshrunk$Q$rho_n))
})


test_that("full-library logCPM preserves matrix identifiers", {
  counts <- Matrix::Matrix(
    matrix(
      c(10, 5, 0, 5),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("G1", "G2"), c("m1", "m2"))
    ),
    sparse = TRUE
  )
  observed <- .rc_metacell_logcpm(
    counts,
    library_size = c(m1 = 100, m2 = 200)
  )
  expect_identical(dimnames(observed), dimnames(counts))
})


test_that("differential helper uses the primary penalty matrix", {
  unit_ids <- paste0("u", 1:4)
  row_id <- "reaction=R1::direction=forward::medium=base"
  result <- list(
    score = matrix(
      c(0.1, 0.2, 0.9, 1.0),
      nrow = 1,
      dimnames = list(row_id, unit_ids)
    ),
    penalty = matrix(
      c(4, 5, 1, 2),
      nrow = 1,
      dimnames = list(row_id, unit_ids)
    ),
    unit_meta = data.frame(
      unit_id = unit_ids,
      sample_id = paste0("s", 1:4),
      condition = c("A", "A", "B", "B"),
      cell_type = "T",
      stringsAsFactors = FALSE
    )
  )

  output <- rc_test_microcompass_differential(
    result,
    method = "lm",
    min_samples_per_group = 2,
    preferred_min_samples_per_group = 2,
    strict_replicate_design = TRUE
  )

  expect_equal(unique(output$analysis_metric), "penalty")
  expect_equal(
    unique(output$metric_direction),
    "positive_effect_means_higher_penalty_and_weaker_support"
  )
  expect_lt(output$effect_size, 0)
})
