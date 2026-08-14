test_that("Pando and SuperCell compatibility is API-based rather than version-based", {
  route <- paste(
    deparse(body(RegCompassR:::.rc_route_pando_infer_args)),
    collapse = "\n"
  )
  expect_match(route, "getNamespaceExports", fixed = TRUE)
  expect_false(grepl("packageVersion", route, fixed = TRUE))

  description <- readLines(testthat::test_path("..", "..", "DESCRIPTION"))
  expect_false(any(grepl("Pando \\(>=", description)))
  expect_false(any(grepl("SuperCell \\(>=", description)))
  remotes <- paste(description, collapse = "\n")
  expect_match(remotes, "1667857557/Pando_regcompass", fixed = TRUE)
  expect_match(remotes, "1667857557/SuperCell_Seurat_V4", fixed = TRUE)
  expect_false(grepl("1667857557/Pando_regcompass@", remotes, fixed = TRUE))
  expect_false(grepl("1667857557/SuperCell_Seurat_V4@", remotes, fixed = TRUE))
})

test_that("condition penalty entry uses Pando activity, fit status and target R2", {
  estimate <- c(0.05, 0.0499, -0.06, 0.06, NA_real_, -0.05)
  padj <- c(0.049, 0.001, 0.049, 0.05, 0.001, 0.001)
  estimable <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
  statistically_supported <- estimable & is.finite(estimate) &
    is.finite(padj) & padj < 0.05
  coefficient <- data.frame(
    estimate = estimate,
    padj = padj,
    estimable = estimable,
    statistically_supported = statistically_supported,
    global_support = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    local_support = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    active = statistically_supported,
    significant = statistically_supported,
    penalty_effect = ifelse(statistically_supported, estimate, 0),
    fit_status = c("ok", "ok", "ok", "ok", "ok", "ok"),
    rsq = c(0.8, 0.049, 0.2, 0.9, 0.9, 0.9),
    stringsAsFactors = FALSE
  )
  expect_identical(
    RegCompassR:::.rc_condition_penalty_gate(
      coefficient, padj_threshold = 0.05, target_rsq_threshold = 0.05
    ),
    c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE)
  )
})

test_that("standard Pando active-edge BH and target R2 thresholds are configurable", {
  edges <- data.frame(
    estimate = c(0.05, 0.049, -0.06, 0.06, -0.05),
    padj = c(0.01, 0.03, 0.05, 0.08, 0.049),
    rsq = c(0.8, 0.2, 0.9, 0.9, 0.9),
    estimable = c(TRUE, TRUE, TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  selected_005 <- RegCompassR:::.rc_filter_standard_pando_edges(
    edges, padj_threshold = 0.05, target_rsq_threshold = 0.05
  )
  selected_002 <- RegCompassR:::.rc_filter_standard_pando_edges(
    edges, padj_threshold = 0.02, target_rsq_threshold = 0.05
  )
  selected_rsq <- RegCompassR:::.rc_filter_standard_pando_edges(
    edges, padj_threshold = 0.05, target_rsq_threshold = 0.5
  )
  expect_equal(selected_005$estimate, c(0.05, 0.049))
  expect_equal(selected_002$estimate, 0.05)
  expect_equal(selected_rsq$estimate, 0.05)
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
