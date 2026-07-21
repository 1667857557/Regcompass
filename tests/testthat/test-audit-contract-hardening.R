test_that("v1.7.0 workflow fixes the canonical GPR architecture", {
  workflow_text <- paste(
    deparse(body(.rc_run_regcompass_uncorrected_metadata)),
    collapse = "\n"
  )
  layer1_text <- paste(
    deparse(body(.rc_build_condition_pooled_layer1)),
    collapse = "\n"
  )

  expect_match(workflow_text, "samples_mixed_within_condition", fixed = TRUE)
  expect_match(layer1_text, 'promiscuity_mode = "none"', fixed = TRUE)
  expect_match(layer1_text, 'and_method = "boltzmann"', fixed = TRUE)
  expect_match(layer1_text, 'or_method = "sum"', fixed = TRUE)
  expect_match(layer1_text, "gene_support_multiome", fixed = TRUE)
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

test_that("condition comparison uses the primary penalty under shared vmax", {
  row_id <- "reaction=R1::direction=forward::medium=base"
  result <- list(
    penalty = matrix(
      c(4, 5, 1, 2),
      nrow = 1,
      dimnames = list(row_id, paste0("u", 1:4))
    ),
    vmax = matrix(
      1,
      nrow = 1,
      ncol = 4,
      dimnames = list(row_id, paste0("u", 1:4))
    ),
    params = list(omega = 1),
    unit_meta = data.frame(
      unit_id = paste0("u", 1:4),
      condition = c("A", "A", "B", "B"),
      cell_type = "T",
      stringsAsFactors = FALSE
    )
  )
  out <- .rc_condition_penalty_comparison(result)
  expect_equal(nrow(out$contrast), 1)
  expect_equal(out$contrast$higher_supported_condition, "B")
  expect_gt(out$contrast$delta_support_b_minus_a, 0)
})
