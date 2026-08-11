test_that("Pando compatibility is API-based rather than version-based", {
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_require_condition_pando_version)),
    collapse = "\n"
  )
  expect_match(implementation, "getNamespaceExports", fixed = TRUE)
  expect_false(grepl("packageVersion", implementation, fixed = TRUE))
  repository_check <- paste(
    deparse(body(RegCompassR:::.rc_validate_pando_repository)),
    collapse = "\n"
  )
  expect_match(repository_check, "required_exported_api_without_version_floor",
               fixed = TRUE)
  expect_false(grepl("package_version", repository_check, fixed = TRUE))
  description <- readLines(testthat::test_path("..", "..", "DESCRIPTION"))
  expect_false(any(grepl("Pando \\(>=", description)))
})

test_that("condition penalty entry uses fit validity and configured BH gate", {
  coefficient <- data.frame(
    estimate = c(0.05, 0.0499, -0.06, 0.06, NA_real_, -0.05),
    padj = c(0.049, 0.001, 0.049, 0.05, 0.001, 0.001),
    estimable = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
    fit_status = c("ok", "ok", "ok", "ok", "ok", "ok"),
    stringsAsFactors = FALSE
  )
  expect_identical(
    RegCompassR:::.rc_condition_penalty_gate(
      coefficient, padj_threshold = 0.05
    ),
    c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)
  )
})

test_that("standard Pando active-edge BH threshold is configurable", {
  edges <- data.frame(
    estimate = c(0.05, 0.049, -0.06, 0.06, -0.05),
    padj = c(0.01, 0.03, 0.05, 0.08, 0.049),
    estimable = c(TRUE, TRUE, TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  selected_005 <- RegCompassR:::.rc_filter_standard_pando_edges(
    edges, padj_threshold = 0.05
  )
  selected_002 <- RegCompassR:::.rc_filter_standard_pando_edges(
    edges, padj_threshold = 0.02
  )
  expect_equal(selected_005$estimate, c(0.05, 0.049))
  expect_equal(selected_002$estimate, 0.05)
  expect_identical(attr(selected_005, "edge_filter")$padj, "< 0.05")
  expect_identical(attr(selected_002, "edge_filter")$padj, "< 0.02")
})

test_that("Gurobi is pinned to one thread inside every worker", {
  implementation <- paste(deparse(body(rc_solve_lp)), collapse = "\n")
  expect_match(implementation, "Threads = 1", fixed = TRUE)
  expect_match(implementation, ".rc_solve_lp_gurobi_base", fixed = TRUE)
})

test_that("FASTCORE and directional LPs expose reaction-granular tasks", {
  fastcore <- paste(
    deparse(body(RegCompassR:::.rc_build_celltype_medium_union_gem_cache)),
    collapse = "\n"
  )
  feasibility <- paste(
    deparse(body(RegCompassR:::.rc_directional_feasibility)),
    collapse = "\n"
  )
  vmax <- paste(
    deparse(body(RegCompassR:::.rc_build_microcompass_vmax_cache)),
    collapse = "\n"
  )
  engine <- paste(
    deparse(body(RegCompassR:::.rc_run_celltype_microcompass_engine)),
    collapse = "\n"
  )
  expect_match(fastcore, "expand.grid", fixed = TRUE)
  expect_match(fastcore, "rc_parallel_lapply", fixed = TRUE)
  expect_match(fastcore, "cell_type_x_medium", fixed = TRUE)
  expect_match(feasibility, "rc_parallel_lapply", fixed = TRUE)
  expect_match(feasibility, ".rc_layer2_task_bpparam", fixed = TRUE)
  expect_match(vmax, "as.list(row_ids)", fixed = TRUE)
  expect_match(vmax, "vmax_reaction_cache", fixed = TRUE)
  expect_match(engine,
    "directional_reaction_by_matching_metacells_step2", fixed = TRUE)
  expect_match(engine, ".rc_atomic_save_rds", fixed = TRUE)
  expect_match(engine, "gc(verbose = FALSE, full = TRUE)", fixed = TRUE)
})