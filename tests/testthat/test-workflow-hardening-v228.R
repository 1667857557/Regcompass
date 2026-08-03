test_that("Stage 1 exposes the fixed min_cells and fragment policy", {
  expect_identical(RegCompassR:::.rc_stage1_min_cells_fixed, 300L)
  expect_identical(formals(rc_regcompass_step_grn)$fragment_files, FALSE)
  body_text <- paste(c(
    deparse(body(rc_regcompass_step_grn)),
    deparse(body(RegCompassR:::.rc_build_stage_analysis_cell_set))
  ), collapse = "\n")
  expect_match(body_text, ".rc_resolve_stage1_min_cells_contract", fixed = TRUE)
  expect_match(body_text, ".rc_filter_stage1_groups_by_min_cells", fixed = TRUE)
  expect_match(body_text, ".rc_drop_zero_count_atac_features", fixed = TRUE)
  expect_match(body_text, ".rc_clear_signac_fragments", fixed = TRUE)
  expect_match(body_text, "stage1_exact_cell_ids_v1", fixed = TRUE)
  filter_text <- paste(
    deparse(body(RegCompassR:::.rc_filter_stage1_groups_by_min_cells)),
    collapse = "\n"
  )
  expect_match(filter_text, "condition_x_cell_type", fixed = TRUE)
  expect_match(filter_text, "retained_stratum", fixed = TRUE)
  expect_match(filter_text, "retained_condition_count >= 2L", fixed = TRUE)
})

test_that("standard Pando uses strict adjusted-P and coefficient filters", {
  expect_identical(RegCompassR:::.rc_standard_pando_padj_fixed, 0.05)
  expect_identical(RegCompassR:::.rc_standard_pando_min_abs_fixed, 0.01)
  body_text <- paste(
    deparse(body(RegCompassR:::.rc_filter_standard_pando_edges)),
    collapse = "\n"
  )
  expect_match(body_text, "padj < .rc_standard_pando_padj_fixed", fixed = TRUE)
  expect_match(body_text, "abs(estimate) > abs_threshold", fixed = TRUE)
  fit_text <- paste(
    deparse(body(RegCompassR:::.rc_fit_standard_pando_by_cell_type)),
    collapse = "\n"
  )
  expect_match(fit_text, "require_padj = TRUE", fixed = TRUE)
  expect_match(fit_text, ".rc_filter_standard_pando_edges", fixed = TRUE)
})

test_that("cell-type TF-IDF uses one triplet reconstruction", {
  body_text <- paste(
    deparse(body(RegCompassR:::.rc_apply_celltype_shared_tfidf)),
    collapse = "\n"
  )
  expect_match(body_text, "TsparseMatrix", fixed = TRUE)
  expect_match(body_text, "Matrix::sparseMatrix", fixed = TRUE)
  expect_match(body_text, "single_triplet_reconstruction", fixed = TRUE)
  expect_false(grepl("output[keep, ] <-", body_text, fixed = TRUE))
})

test_that("every parallel task sees a one-thread solver environment", {
  old <- Sys.getenv("HIGHS_THREADS", unset = NA_character_)
  on.exit({
    if (is.na(old)) Sys.unsetenv("HIGHS_THREADS") else
      Sys.setenv(HIGHS_THREADS = old)
  }, add = TRUE)
  observed <- rc_parallel_lapply(
    list(1L),
    function(x) Sys.getenv("HIGHS_THREADS", unset = "missing"),
    BPPARAM = FALSE
  )
  expect_identical(observed[[1L]], "1")
})

test_that("Layer 2 propagates progress through direct engine calls", {
  body_text <- paste(deparse(body(rc_regcompass_step_layer2)), collapse = "\n")
  expect_match(body_text, ".rc_run_microcompass_monitored", fixed = TRUE)
  expect_match(body_text, ".rc_run_microcompass_engine_monitored", fixed = TRUE)
  expect_match(body_text, ".rc_compact_meta_modules_for_layer2", fixed = TRUE)
  expect_false(grepl("evaluation_environment", body_text, fixed = TRUE))
})

test_that("external condition-module references are checksummed and loadable", {
  payload <- list(
    reaction_membership = data.frame(
      reaction_id = "R1", stringsAsFactors = FALSE
    )
  )
  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file, force = TRUE), add = TRUE)
  saveRDS(payload, file)
  reference <- list(
    schema_version = "regcompass_external_condition_modules_v1",
    embedded = FALSE,
    file = file,
    file_checksum = unname(tools::md5sum(file))
  )
  meta_modules <- list(
    condition_modules = reference,
    condition_modules_ref = reference
  )
  expect_identical(
    RegCompassR:::.rc_load_condition_modules(meta_modules),
    payload
  )
})

test_that("Stage 3 stores a lightweight condition-module reference", {
  body_text <- paste(
    deparse(body(rc_regcompass_step_meta_modules)),
    collapse = "\n"
  )
  expect_match(body_text, "regcompass_external_condition_modules_v1", fixed = TRUE)
  expect_match(body_text, "condition_modules_embedded = FALSE", fixed = TRUE)
  expect_match(body_text, "condition_modules_ref = condition_modules_ref", fixed = TRUE)
})

test_that("workflow code has no late-loaded override chain", {
  root <- testthat::test_path("..", "..")
  expect_false(any(file.exists(file.path(root, "R", c(
    "zzz_workflow_hardening.R",
    "zzzz_min_cells_contract.R",
    "zzzzz_stage_cell_set_contract.R",
    "zzzzzz_strict_condition_penalty.R"
  )))))
  r_files <- list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)
  source <- unlist(lapply(r_files, readLines, warn = FALSE), use.names = FALSE)
  expect_false(any(grepl("\\.rc_original_.*hardening|evaluation_environment", source)))
})