test_that("RegCompass declares the common-dictionary Pando dependency", {
  package_root <- testthat::test_path("..", "..")
  description <- read.dcf(file.path(package_root, "DESCRIPTION"))
  expect_identical(unname(description[1L, "Version"]), "2.3.0")
  expect_identical(
    unname(description[1L, "Config/RegCompass/PandoConditionSchema"]),
    "pando_condition_grn_common_dictionary_v1"
  )
  expect_identical(
    unname(description[1L, "Config/RegCompass/PandoConditionMethod"]),
    "two-stage-exact-edge-union-fixed-dictionary-glm"
  )
  expect_identical(
    unname(description[1L, "Config/RegCompass/PandoConditionEffectFilter"]),
    "BH-adjusted-p-below-0.05"
  )
  remotes <- gsub(
    "[[:space:]]+", " ",
    trimws(unname(description[1L, "Remotes"]))
  )
  expect_match(
    remotes,
    "1667857557/Pando_regcompass@agent/common-dictionary-condition-grn",
    fixed = TRUE
  )
  collate <- unname(description[1L, "Collate"])
  expect_match(collate, "step_grn_common_dictionary.R", fixed = TRUE)
  expect_match(collate, "step_layer1_common_dictionary.R", fixed = TRUE)
  expect_match(
    collate, "step_layer1_common_dictionary_contract.R", fixed = TRUE
  )
  expect_false(grepl("condition_grn_runtime.R", collate, fixed = TRUE))
  expect_false(grepl("condition_grn_runtime_guard.R", collate, fixed = TRUE))
  expect_false(file.exists(file.path(
    package_root, "R", "condition_grn_runtime.R"
  )))
})

test_that("an installed Pando exposes the common-dictionary API", {
  skip_if_not_installed("Pando")
  description <- utils::packageDescription("Pando")
  expect_identical(
    description[["Config/Pando/ConditionGRNSchema"]],
    "pando_condition_grn_common_dictionary_v1"
  )
  expect_identical(
    description[["Config/Pando/ConditionGRNMethod"]],
    "two-stage-exact-edge-union-fixed-dictionary-glm"
  )
  namespace <- asNamespace("Pando")
  required <- c(
    "infer_condition_grn",
    "condition_grn_fit",
    "condition_grn_subgraph",
    "discover_grn_edges",
    "union_grn_edges",
    "fit_grn_from_edges",
    "project_condition_grn_cells",
    "aggregate_condition_grn_projection"
  )
  expect_true(all(vapply(required, exists, logical(1),
                         envir = namespace, inherits = FALSE)))
  retired <- c(
    ".condition_fit_target_engine_cpp",
    ".condition_fit_multitask_path_cpp",
    ".condition_refit_path_cpp",
    ".condition_product_matrix_cpp"
  )
  expect_false(any(vapply(retired, exists, logical(1),
                          envir = namespace, inherits = FALSE)))
})
