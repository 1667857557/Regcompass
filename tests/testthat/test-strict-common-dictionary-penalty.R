.current_edge_fixture <- function(beta = matrix(
    c(0.4, 0.1, 0.02, 0.01), nrow = 2L,
    dimnames = list(c("A", "B"), c("E1", "E2")))) {
  conditions <- rownames(beta)
  edges <- colnames(beta)
  variance <- 0.01
  edge_stat <- colSums(beta^2 / variance)
  edge_p <- stats::pchisq(edge_stat, df = nrow(beta), lower.tail = FALSE)
  edge_q <- stats::p.adjust(edge_p, method = "BH")
  rows <- lapply(seq_along(edges), function(j) {
    b <- beta[, j]
    stat <- b / sqrt(variance)
    supported <- edge_q[[j]] < 0.05
    data.frame(
      edge_id = edges[[j]], target = "G", condition = conditions,
      estimate = 0.8 + 0.1 * j, penalty_effect = 0.8 + 0.1 * j,
      inference_estimate = b,
      inference_se = sqrt(variance),
      inference_variance = variance,
      inference_statistic = stat,
      condition_pval = 2 * stats::pnorm(-abs(stat)),
      condition_inference_estimable = TRUE,
      edge_df = as.integer(nrow(beta)),
      edge_statistic = edge_stat[[j]],
      edge_pval = edge_p[[j]], edge_padj = edge_q[[j]],
      edge_inference_estimable = TRUE,
      edge_inference_test = "independent_condition_wald_chisq",
      pval = edge_p[[j]], padj = edge_q[[j]],
      bh_scope = "exact_edge_whole_cell_type_network_BH",
      bh_family_size = length(edges),
      all_conditions_fit_valid = TRUE,
      edge_supported = supported,
      statistically_supported = supported,
      significant = supported,
      pando_estimation_active = TRUE, active = TRUE,
      active_in_regcompass = supported,
      fit_status = "ok",
      penalty_family = "information_scaled_sparse_deviation",
      penalty_value = 0.25, solver_status = "ok",
      kkt_residual = 1e-10, iterations = 20L,
      contrast_identifiable = c(FALSE, FALSE),
      shared_by_boundary = c(TRUE, TRUE),
      fused_by_penalty = c(FALSE, FALSE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

test_that("conditional penalty validates exact-edge omnibus and one network BH", {
  coefficient <- .current_edge_fixture()
  expect_invisible(
    RegCompassR:::.rc_validate_pando_active_condition_edges(
      coefficient, padj_threshold = 0.05
    )
  )
  expect_identical(
    unique(coefficient$bh_family_size),
    2L
  )
  expect_true(all(
    coefficient$bh_scope == "exact_edge_whole_cell_type_network_BH"
  ))

  first <- coefficient$edge_id == "E1"
  second <- coefficient$edge_id == "E2"
  expect_true(all(coefficient$active_in_regcompass[first]))
  expect_false(any(coefficient$active_in_regcompass[second]))
})

test_that("condition-local power cannot create condition-specific topology", {
  coefficient <- .current_edge_fixture()
  first <- which(coefficient$edge_id == "E1")
  expect_equal(length(unique(coefficient$condition_pval[first])), 2L)
  expect_true(all(coefficient$active_in_regcompass[first]))

  broken <- coefficient
  broken$active_in_regcompass[first[[2L]]] <- FALSE
  expect_error(
    RegCompassR:::.rc_validate_pando_active_condition_edges(
      broken, padj_threshold = 0.05
    ),
    "topology"
  )
})

test_that("boundary sharing is production metadata rather than significance", {
  coefficient <- .current_edge_fixture()
  gate <- RegCompassR:::.rc_condition_penalty_gate(
    coefficient, padj_threshold = 0.05,
    target_rsq_threshold = 0.95
  )
  expect_identical(gate, coefficient$active_in_regcompass)
  expect_true(all(coefficient$shared_by_boundary))
  expect_false(any(coefficient$contrast_identifiable))
})

test_that("failed production invalidates an otherwise supported common edge", {
  coefficient <- .current_edge_fixture()
  first <- which(coefficient$edge_id == "E1")
  failed <- coefficient
  failed$fit_status[first[[2L]]] <- "failed"
  failed$all_conditions_fit_valid[first] <- FALSE
  failed$edge_supported[first] <- FALSE
  failed$statistically_supported[first] <- FALSE
  failed$significant[first] <- FALSE
  failed$active_in_regcompass[first] <- FALSE

  gate <- RegCompassR:::.rc_condition_penalty_gate(
    failed, padj_threshold = 0.05,
    target_rsq_threshold = 0.05
  )
  expect_false(any(gate[first]))
})

test_that("Layer 1 retains canonical RegCompass TF-ATAC exposure mixture", {
  helper <- paste(
    deparse(body(RegCompassR:::.rc_pando_projection_from_group_means)),
    collapse = "\n"
  )
  projection <- paste(
    deparse(body(RegCompassR:::.rc_condition_pando_projection)),
    collapse = "\n"
  )
  expect_match(helper, "paired_product <- tf_block * peak_block", fixed = TRUE)
  expect_match(
    helper,
    "product_of_means_weight * tf_mean * peak_mean",
    fixed = TRUE
  )
  expect_match(projection, ".rc_condition_penalty_gate", fixed = TRUE)
  expect_identical(
    RegCompassR:::.RC_PANDO_PROJECTION_PRODUCT_OF_MEANS_WEIGHT,
    0.25
  )
})
