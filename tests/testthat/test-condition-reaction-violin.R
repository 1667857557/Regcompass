make_condition_violin_fixture <- function() {
  conditions <- rep(c("Control", "Treatment"), each = 6L)
  units <- paste0("v", seq_along(conditions))
  row_id <- "reaction=R_violin::direction=forward::medium=base"
  penalty <- matrix(
    c(12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1),
    nrow = 1L,
    dimnames = list(row_id, units)
  )
  list(
    penalty = penalty,
    vmax = matrix(100, 1L, length(units), dimnames = dimnames(penalty)),
    feasible = matrix(TRUE, 1L, length(units), dimnames = dimnames(penalty)),
    unit_meta = data.frame(
      pool_id = units,
      condition = conditions,
      cell_type = "T_cell",
      stringsAsFactors = FALSE
    ),
    params = list(omega = 0.95, unit = "metacell")
  )
}

test_that("selected reaction plots retain violin geometry", {
  skip_if_not_installed("ggplot2")
  plot <- rc_plot_condition_reaction(
    make_condition_violin_fixture(),
    reaction_id = "R_violin",
    cell_type = "T_cell",
    target_direction = "forward",
    medium_scenario = "base",
    condition_col = "condition",
    celltype_col = "cell_type",
    conditions = c("Control", "Treatment"),
    plot_type = "violin_boxplot",
    show_nonsignificant = TRUE
  )

  expect_s3_class(plot, "ggplot")
  expect_identical(attr(plot, "plot_type"), "violin_boxplot")
  geom_classes <- vapply(plot$layers, function(layer) {
    class(layer$geom)[[1L]]
  }, character(1))
  expect_true("GeomViolin" %in% geom_classes)
  expect_true("GeomBoxplot" %in% geom_classes)
})

test_that("violin-only reaction plots remove the boxplot layer", {
  skip_if_not_installed("ggplot2")
  plot <- rc_plot_condition_reaction(
    make_condition_violin_fixture(),
    reaction_id = "R_violin",
    cell_type = "T_cell",
    target_direction = "forward",
    medium_scenario = "base",
    condition_col = "condition",
    celltype_col = "cell_type",
    conditions = c("Control", "Treatment"),
    plot_type = "violin"
  )

  geom_classes <- vapply(plot$layers, function(layer) {
    class(layer$geom)[[1L]]
  }, character(1))
  expect_true("GeomViolin" %in% geom_classes)
  expect_false("GeomBoxplot" %in% geom_classes)
})

test_that("invalid violin controls fail before plotting", {
  skip_if_not_installed("ggplot2")
  expect_error(
    rc_plot_condition_reaction(
      make_condition_violin_fixture(),
      reaction_id = "R_violin",
      cell_type = "T_cell",
      target_direction = "forward",
      medium_scenario = "base",
      violin_alpha = 2
    ),
    "Violin geometry controls"
  )
})
