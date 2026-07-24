test_that("GRN and metacell groups require bidirectional coverage", {
  grn <- list(sample_status = data.frame(
    condition = c("A", "B"),
    cell_type = c("T", "T"),
    status = c("ok", "ok"),
    n_cells = c(100L, 120L),
    n_significant_edges = c(10L, 12L),
    stringsAsFactors = FALSE
  ))
  metacells <- data.frame(
    metacell_id = c("A1", "A2", "B1"),
    condition = c("A", "A", "B"),
    cell_type = "T",
    dominant_celltype_fraction = c(1, 0.8, 1),
    mixed_celltype_metacell = c(FALSE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  coverage <- .rc_validate_grn_metacell_group_coverage(
    grn, metacells, "condition", "cell_type"
  )
  expect_true(all(coverage$coverage_complete))
  expect_equal(coverage$n_metacells[coverage$condition == "A"], 2L)
  expect_equal(
    coverage$n_mixed_celltype_metacells[coverage$condition == "A"],
    1L
  )

  expect_error(
    .rc_validate_grn_metacell_group_coverage(
      grn,
      metacells[metacells$condition == "A", , drop = FALSE],
      "condition", "cell_type"
    ),
    "do not align"
  )
})

test_that("condition-only metacells reject tied dominant cell types", {
  skip_if_not_installed("SeuratObject")
  counts <- Matrix::Matrix(
    matrix(
      seq_len(8), nrow = 2,
      dimnames = list(c("g1", "g2"), paste0("c", 1:4))
    ),
    sparse = TRUE
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object$cell_type <- c("T", "B", "T", "B")
  pooled <- list(
    membership = data.frame(
      cell_id = paste0("c", 1:4),
      metacell_id = c("mc1", "mc1", "mc2", "mc2"),
      stringsAsFactors = FALSE
    ),
    metacell_meta = data.frame(
      metacell_id = c("mc1", "mc2"),
      stringsAsFactors = FALSE
    )
  )

  expect_error(
    .rc_assign_metacell_dominant_celltype(pooled, object, "cell_type"),
    "tied dominant cell types"
  )
})

test_that("metacell stage persists post hoc metadata contracts", {
  text <- paste(deparse(body(rc_regcompass_step_metacells)), collapse = "\n")
  required <- c(
    "metacell_metadata.tsv.gz",
    "metacell_membership.tsv.gz",
    "metacell_celltype_composition.tsv.gz",
    "metacell_celltype_summary.tsv.gz",
    "merged_metacell_object.rds",
    "step_metacells.rds"
  )
  expect_true(all(vapply(required, grepl, logical(1), x = text, fixed = TRUE)))
})

test_that("canonical Layer 1 has no versioned compatibility override", {
  description <- utils::packageDescription("RegCompassR")
  collate <- description$Collate %||% ""
  expect_false(grepl("v170_layer1_parallel.R", collate, fixed = TRUE))
  body_text <- paste(
    deparse(body(.rc_build_condition_pooled_layer1)), collapse = "\n"
  )
  expect_match(
    body_text,
    "regcompass_condition_only_layer1_v1.8.1",
    fixed = TRUE
  )
  expect_match(
    body_text,
    "condition_only_metacell_with_posthoc_celltype",
    fixed = TRUE
  )
})

test_that("final results retain modules and add reaction interpretation", {
  body_text <- paste(deparse(body(rc_regcompass_step_results)), collapse = "\n")
  expect_match(body_text, "condition_grn_meta_modules", fixed = TRUE)
  expect_match(body_text, "global_grn_meta_modules", fixed = TRUE)
  expect_match(body_text, "grn_metacell_group_coverage", fixed = TRUE)
  expect_match(body_text, "reaction_catalog", fixed = TRUE)
  expect_match(body_text, "reaction_evidence", fixed = TRUE)
})
