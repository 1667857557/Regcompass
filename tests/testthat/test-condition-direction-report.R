make_condition_direction_report_fixture <- function() {
  conditions <- rep(c("control", "JQ1", "MS177"), each = 6L)
  units <- paste0("u", seq_along(conditions))
  row_ids <- c(
    "reaction=R_symmetric::direction=forward::medium=base",
    "reaction=R_symmetric::direction=reverse::medium=base",
    "reaction=R_asymmetric::direction=forward::medium=base",
    "reaction=R_asymmetric::direction=reverse::medium=base",
    "reaction=R_forward_only::direction=forward::medium=base"
  )
  shifted_penalty <- c(
    10, 11, 12, 13, 14, 15,
    5, 6, 7, 8, 9, 10,
    1, 2, 3, 4, 5, 6
  )
  penalty <- rbind(
    R_symmetric_forward = shifted_penalty,
    R_symmetric_reverse = shifted_penalty,
    R_asymmetric_forward = shifted_penalty,
    R_asymmetric_reverse = rep(8, length(units)),
    R_forward_only = shifted_penalty
  )
  rownames(penalty) <- row_ids
  colnames(penalty) <- units
  list(
    penalty = penalty,
    vmax = matrix(
      100,
      nrow = nrow(penalty),
      ncol = ncol(penalty),
      dimnames = dimnames(penalty)
    ),
    feasible = matrix(
      TRUE,
      nrow = nrow(penalty),
      ncol = ncol(penalty),
      dimnames = dimnames(penalty)
    ),
    unit_meta = data.frame(
      pool_id = units,
      condition = conditions,
      cell_type = "epithelial_like",
      stringsAsFactors = FALSE
    ),
    params = list(omega = 0.95, unit = "metacell")
  )
}

test_that("direction report preserves directional results and adds summaries", {
  report <- rc_report_condition_directions(
    make_condition_direction_report_fixture(),
    condition_col = "condition",
    celltype_col = "cell_type",
    min_units = 5L,
    include_unit_metrics = TRUE,
    source_label = "fixture"
  )

  expect_s3_class(report, "regcompass_condition_direction_report")
  expect_equal(nrow(report$directional_pairwise), 15L)
  expect_equal(nrow(report$directional_omnibus), 5L)
  expect_equal(nrow(report$reaction_pairwise), 15L)
  expect_equal(nrow(report$reaction_omnibus), 5L)
  expect_equal(nrow(report$direction_diagnostics), 9L)
  expect_equal(nrow(report$unit_metrics), 54L)
  expect_true(all(report$directional_pairwise$target_direction %in%
    c("forward", "reverse")))
  expect_true(all(report$reaction_pairwise$report_metric %in%
    c("any_direction_support", "directional_balance")))
  expect_match(report$reporting_policy, "not net flux")
})

test_that("numerically identical directions are reported without double counting", {
  report <- rc_report_condition_directions(
    make_condition_direction_report_fixture(),
    condition_col = "condition",
    celltype_col = "cell_type",
    min_units = 5L,
    direction_tolerance = 1e-12
  )

  symmetric <- subset(
    report$direction_diagnostics,
    reaction_id == "R_symmetric"
  )
  expect_true(all(symmetric$directionally_indistinguishable))
  expect_true(all(
    symmetric$direction_pair_status == "bidirectional_indistinguishable"
  ))
  expect_true(all(symmetric$preferred_direction == "indistinguishable"))
  expect_equal(symmetric$median_directional_balance, rep(0, 3L))

  balance <- subset(
    report$reaction_pairwise,
    reaction_id == "R_symmetric" &
      report_metric == "directional_balance"
  )
  expect_equal(balance$delta_median_metric_b_minus_a, rep(0, 3L))
  expect_equal(balance$p_value, rep(1, 3L))

  any_direction <- subset(
    report$reaction_pairwise,
    reaction_id == "R_symmetric" &
      report_metric == "any_direction_support" &
      condition_a == "control" & condition_b == "MS177"
  )
  expect_equal(nrow(any_direction), 1L)
  expect_gt(any_direction$delta_median_metric_b_minus_a, 0)
})

test_that("directional balance reports support shifts but not net flux", {
  report <- rc_report_condition_directions(
    make_condition_direction_report_fixture(),
    condition_col = "condition",
    celltype_col = "cell_type",
    min_units = 5L
  )

  asymmetric <- subset(
    report$direction_diagnostics,
    reaction_id == "R_asymmetric"
  )
  expect_false(any(asymmetric$directionally_indistinguishable))

  shift <- subset(
    report$reaction_pairwise,
    reaction_id == "R_asymmetric" &
      report_metric == "directional_balance" &
      condition_a == "control" & condition_b == "MS177"
  )
  expect_equal(nrow(shift), 1L)
  expect_gt(shift$delta_median_metric_b_minus_a, 0)
  expect_equal(shift$direction_shift_b_minus_a, "toward_forward")
  expect_match(shift$metric_semantics, "not net flux")

  forward_only <- subset(
    report$direction_diagnostics,
    reaction_id == "R_forward_only"
  )
  expect_true(all(forward_only$direction_pair_status == "forward_only"))
  expect_false(any(
    report$reaction_pairwise$reaction_id == "R_forward_only" &
      report$reaction_pairwise$report_metric == "directional_balance"
  ))
})

test_that("direction report exports all report layers", {
  outdir <- tempfile("condition_direction_report_")
  report <- rc_report_condition_directions(
    make_condition_direction_report_fixture(),
    condition_col = "condition",
    celltype_col = "cell_type",
    min_units = 5L,
    include_unit_metrics = TRUE,
    outdir = outdir
  )

  expect_true(nrow(report$reaction_pairwise) > 0L)
  expected <- c(
    "condition_directional_pairwise.tsv.gz",
    "condition_directional_omnibus.tsv.gz",
    "condition_reaction_level_pairwise.tsv.gz",
    "condition_reaction_level_omnibus.tsv.gz",
    "condition_direction_diagnostics.tsv.gz",
    "condition_direction_unit_metrics.tsv.gz",
    "condition_direction_report.rds"
  )
  expect_true(all(file.exists(file.path(outdir, expected))))
})
