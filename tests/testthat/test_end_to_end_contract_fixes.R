test_that("HiGHS numeric status 7 with valid fields is optimal", {
  result <- list(
    status = 7L,
    model_status = "Optimal",
    objective_value = 1,
    primal_solution = c(1, 0)
  )
  expect_equal(
    .rc_highs_status(result, n_variables = 2),
    "optimal"
  )
  expect_equal(
    .rc_highs_status(
      list(
        status = "optimal",
        objective_value = 1,
        primal_solution = c(1, 0)
      ),
      n_variables = 2
    ),
    "optimal"
  )
})

test_that("infeasible COMPASS scores are NA rather than zero", {
  penalty <- matrix(
    c(1, 2),
    nrow = 1,
    dimnames = list("R1", c("u1", "u2"))
  )
  feasible <- matrix(
    c(TRUE, FALSE),
    nrow = 1,
    dimnames = dimnames(penalty)
  )
  score <- rc_compass_score_from_penalty(
    penalty,
    feasible
  )
  expect_true(is.na(score["R1", "u2"]))
})

test_that("single-sample output is descriptive-only when requested", {
  result <- list(
    score = matrix(
      c(0.2, 0.8),
      nrow = 1,
      dimnames = list(
        "R1::forward::base",
        c("u1", "u2")
      )
    ),
    unit_meta = data.frame(
      unit_id = c("u1", "u2"),
      sample_id = "s1",
      condition = "ctrl",
      cell_type = "T",
      n_cells = c(10, 20),
      stringsAsFactors = FALSE
    )
  )
  expect_error(
    rc_test_microcompass_differential(
      result,
      strict_replicate_design = TRUE
    ),
    "biological samples"
  )
  output <- rc_test_microcompass_differential(
    result,
    strict_replicate_design = FALSE
  )
  expect_true(all(is.na(output$p_value)))
  expect_true(all(is.na(output$FDR)))
  expect_equal(unique(output$medium_scenario), "base")
  expect_equal(unique(output$model_status), "descriptive_only")
  expect_equal(unique(output$n_biological_samples), 1L)
})

test_that("strict replicate design is enforced within cell type", {
  unit_meta <- rbind(
    data.frame(
      unit_id = paste0("T", 1:6),
      sample_id = paste0("s", 1:6),
      condition = rep(c("A", "B"), each = 3),
      cell_type = "T",
      stringsAsFactors = FALSE
    ),
    data.frame(
      unit_id = paste0("B", 1:4),
      sample_id = paste0("s", 7:10),
      condition = c("A", "A", "A", "B"),
      cell_type = "B",
      stringsAsFactors = FALSE
    )
  )
  result <- list(
    score = matrix(
      seq_len(nrow(unit_meta)),
      nrow = 1,
      dimnames = list(
        "R1::forward::base",
        unit_meta$unit_id
      )
    ),
    unit_meta = unit_meta
  )
  expect_error(
    rc_test_microcompass_differential(
      result,
      method = "lm",
      min_samples_per_group = 2,
      strict_replicate_design = TRUE
    ),
    "within cell type"
  )
})

test_that("one shared fragment file maps to every sample", {
  file <- tempfile(fileext = ".tsv.gz")
  writeLines("chr1\t1\t2\tcell-1\t1", file)
  manifest <- .rc_normalize_fragment_manifest(
    file,
    sample_ids = c("sample1", "sample2"),
    atac_assay = "ATAC"
  )
  expect_equal(
    manifest$sample_id,
    c("sample1", "sample2")
  )
  expect_equal(unique(manifest$fragment_file), file)
})

test_that("multi-class condition lm uses an omnibus test", {
  result <- list(
    score = matrix(
      c(1, 2, 3, 1.2, 2.2, 3.2),
      nrow = 1,
      dimnames = list(
        "R1::forward::base",
        paste0("u", 1:6)
      )
    ),
    unit_meta = data.frame(
      unit_id = paste0("u", 1:6),
      sample_id = paste0("s", 1:6),
      condition = rep(c("A", "B", "C"), each = 2),
      cell_type = "T",
      stringsAsFactors = FALSE
    )
  )
  output <- rc_test_microcompass_differential(
    result,
    method = "lm",
    min_samples_per_group = 2,
    strict_replicate_design = TRUE,
    test_type = "omnibus"
  )
  expect_equal(output$contrast, "condition_omnibus")
  expect_true(is.finite(output$p_value))
})
