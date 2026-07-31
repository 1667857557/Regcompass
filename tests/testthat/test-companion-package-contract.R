test_that("companion remotes are immutable audited commits", {
  description <- read.dcf("DESCRIPTION")
  remotes <- description[1L, "Remotes"]
  expect_match(
    remotes,
    "SuperCell_Seurat_V4@1319730b15857d88f169a9006f6cf4e88993c0c7",
    fixed = TRUE
  )
  expect_match(
    remotes,
    "Pando_regcompass@9972d0a576ecc55c8f4898cb85d4fd90f2b30b3b",
    fixed = TRUE
  )
  expect_false(grepl("@Supercell2", remotes, fixed = TRUE))
  expect_false(grepl("@agent/", remotes, fixed = TRUE))
})

test_that("installed Pando retains optimized condition-GRN dispatch", {
  skip_if_not_installed("Pando")
  expect_gte(packageVersion("Pando"), package_version("1.5.0"))
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

test_that("installed SuperCell retains grouped graph semantics", {
  skip_if_not_installed("SuperCell")
  expect_true(
    "SCimplify_by_graph_group_from_embedding" %in%
      getNamespaceExports("SuperCell")
  )
  formals <- names(formals(
    SuperCell::SCimplify_by_graph_group_from_embedding
  ))
  expect_true(all(c(
    "X", "cell.graph.group", "cell.split.condition",
    "gamma", "k.knn", "n.pc"
  ) %in% formals))
})
