test_that("companion remotes install the latest repository heads", {
  description_path <- testthat::test_path("..", "..", "DESCRIPTION")
  description <- read.dcf(description_path)
  remotes <- description[1L, "Remotes"]
  expect_match(
    remotes,
    "1667857557/SuperCell_Seurat_V4",
    fixed = TRUE
  )
  expect_false(grepl("SuperCell_Seurat_V4@", remotes, fixed = TRUE))
  expect_match(
    remotes,
    "1667857557/Pando_regcompass",
    fixed = TRUE
  )
  expect_false(grepl("Pando_regcompass@", remotes, fixed = TRUE))
})

test_that("installed Pando retains optimized condition-GRN dispatch", {
  skip_if_not_installed("Pando")
  expect_true(is.function(Pando::infer_condition_grn))
  expect_true(is.function(Pando::condition_grn_fit))
  expect_true(is.function(Pando::project_condition_grn_primary_cells))
  expect_true(all(
    c("backend", "verified_estimability_mask") %in%
      names(formals(Pando:::.condition_fit_multitask_path))
  ))
  expect_true(all(
    c("cache", "estimability_verified", "solver") %in%
      names(formals(Pando:::.condition_refit_shared_baseline))
  ))
})

test_that("installed SuperCell exposes the canonical grouped WNN API", {
  skip_if_not_installed("SuperCell")
  expect_true("SCimplify_by_graph_group" %in% getNamespaceExports("SuperCell"))
  expect_false(
    "SCimplify_by_graph_group_from_embedding" %in%
      getNamespaceExports("SuperCell")
  )
  api_formals <- names(formals(SuperCell::SCimplify_by_graph_group))
  expect_true(all(c(
    "seurat", "cell.graph.group", "cell.split.condition",
    "assay", "reduction", "dims", "gamma", "k.knn"
  ) %in% api_formals))
  expect_identical(eval(formals(SuperCell::SCimplify_for_Seurat)$gamma), 30)
  expect_identical(eval(formals(SuperCell::SCimplify_by_graph_group)$gamma), 30)
})
