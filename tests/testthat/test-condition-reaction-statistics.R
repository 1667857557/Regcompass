make_condition_statistics_fixture <- function() {
  conditions <- rep(c("control", "JQ1", "MS177"), each = 6L)
  units <- paste0("u", seq_along(conditions))
  row_ids <- c(
    "reaction=R_shift::direction=forward::medium=base",
    "reaction=R_constant::direction=reverse::medium=base"
  )
  shifted_penalty <- c(
    10, 11, 12, 13, 14, 15,
    5, 6, 7, 8, 9, 10,
    1, 2, 3, 4, 5, 6
  )
  penalty <- rbind(
    R_shift = shifted_penalty,
    R_constant = rep(5, length(units))
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

test_that("shared-GEM reaction scores support pairwise condition tests", {
  microcompass <- make_condition_statistics_fixture()
  result <- rc_test_condition_reactions(
    microcompass,
    condition_col = "condition",
    celltype_col = "cell_type",
    min_units = 5L,
    include_scores = TRUE
  )

  expect_s3_class(result, "regcompass_condition_statistics")
  expect_equal(nrow(result$pairwise), 6L)
  expect_equal(nrow(result$omnibus), 2L)
  expect_equal(dim(result$score), c(2L, 18L))

  shifted <- subset(
    result$pairwise,
    reaction_id == "R_shift" &
      condition_a == "control" & condition_b == "MS177"
  )
  expect_equal(nrow(shifted), 1L)
  expect_gt(shifted$delta_median_score_b_minus_a, 0)
  expect_gt(shifted$rank_biserial_b_minus_a, 0)
  expect_equal(shifted$higher_supported_condition, "MS177")
  expect_lt(shifted$p_value, 0.05)
  expect_true(is.finite(shifted$p_adj))
  expect_equal(shifted$inference_level, "metacell_within_dataset")
  expect_true(shifted$descriptive_only)
  expect_false(shifted$biological_replicate_inference)

  constant <- subset(
    result$pairwise,
    reaction_id == "R_constant" &
      condition_a == "control" & condition_b == "JQ1"
  )
  expect_equal(constant$p_value, 1)
  expect_equal(constant$test_status, "constant_equal")
  expect_equal(constant$rank_biserial_b_minus_a, 0)
})

test_that("top-level RegCompass results and filters are accepted", {
  microcompass <- make_condition_statistics_fixture()
  result <- rc_test_condition_reactions(
    list(microcompass = microcompass),
    condition_col = "condition",
    celltype_col = "cell_type",
    conditions = c("control", "MS177"),
    reaction_ids = "R_shift",
    comparisons = list(c("control", "MS177")),
    min_units = 5L
  )
  expect_equal(nrow(result$pairwise), 1L)
  expect_equal(result$pairwise$reaction_id, "R_shift")
  expect_equal(nrow(result$omnibus), 0L)
})

test_that("shared-GEM vmax invariance is enforced", {
  microcompass <- make_condition_statistics_fixture()
  microcompass$vmax[1L, 1L] <- 90
  expect_error(
    rc_test_condition_reactions(
      microcompass,
      condition_col = "condition",
      celltype_col = "cell_type"
    ),
    "vmax differs"
  )
})

test_that("condition statistics can be exported", {
  microcompass <- make_condition_statistics_fixture()
  outdir <- tempfile("condition_statistics_")
  result <- rc_test_condition_reactions(
    microcompass,
    condition_col = "condition",
    celltype_col = "cell_type",
    outdir = outdir
  )
  expect_true(nrow(result$pairwise) > 0L)
  expect_true(file.exists(file.path(
    outdir, "condition_reaction_pairwise.tsv.gz"
  )))
  expect_true(file.exists(file.path(
    outdir, "condition_reaction_omnibus.tsv.gz"
  )))
  expect_true(file.exists(file.path(
    outdir, "condition_reaction_statistics.rds"
  )))
})

test_that("one reaction can be plotted across multiple conditions", {
  skip_if_not_installed("ggplot2")
  microcompass <- make_condition_statistics_fixture()
  plot <- rc_plot_condition_reaction(
    microcompass,
    reaction_id = "R_shift",
    cell_type = "epithelial_like",
    target_direction = "forward",
    medium_scenario = "base",
    condition_col = "condition",
    celltype_col = "cell_type",
    conditions = c("control", "JQ1", "MS177"),
    show_nonsignificant = TRUE
  )

  expect_s3_class(plot, "ggplot")
  plot_data <- attr(plot, "plot_data")
  annotation_data <- attr(plot, "annotation_data")
  statistics <- attr(plot, "condition_statistics")
  expect_equal(nrow(plot_data), 18L)
  expect_equal(
    levels(plot_data$condition),
    c("control", "JQ1", "MS177")
  )
  expect_equal(nrow(annotation_data), 3L)
  expect_true(all(annotation_data$label %in% c("ns", "*", "**", "***", "****")))
  expect_equal(nrow(statistics$pairwise), 6L)
  expect_equal(nrow(statistics$omnibus), 2L)
  expect_true(all(
    statistics$pairwise$p_adj >= statistics$pairwise$p_value - 1e-12,
    na.rm = TRUE
  ))
  built <- ggplot2::ggplot_build(plot)
  expect_gte(length(built$data), 2L)
})