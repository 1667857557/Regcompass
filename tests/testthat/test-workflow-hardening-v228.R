test_that("Stage 1 exposes the fixed prefilter and fragment policy", {
  expect_identical(RegCompassR:::.rc_stage1_min_cells_fixed, 300L)
  expect_identical(formals(rc_regcompass_step_grn)$fragment_files, FALSE)
  body_text <- paste(deparse(body(rc_regcompass_step_grn)), collapse = "\n")
  expect_match(body_text, ".rc_prefilter_stage1_celltypes", fixed = TRUE)
  expect_match(body_text, ".rc_drop_zero_count_atac_features", fixed = TRUE)
  expect_match(body_text, ".rc_clear_signac_fragments", fixed = TRUE)
  expect_match(body_text, "single_cell_grn.rds", fixed = TRUE)
})

test_that("standard Pando uses strict adjusted-P and coefficient filters", {
  expect_identical(RegCompassR:::.rc_standard_pando_padj_fixed, 0.05)
  expect_identical(RegCompassR:::.rc_standard_pando_min_abs_fixed, 0.01)
  body_text <- paste(
    deparse(body(RegCompassR:::.rc_standard_pando_extract_strict)),
    collapse = "\n"
  )
  expect_match(body_text, "padj < .rc_standard_pando_padj_fixed", fixed = TRUE)
  expect_match(body_text, "abs(estimate) > abs_threshold", fixed = TRUE)
  expect_match(body_text, "require_padj <- TRUE", fixed = TRUE)
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

test_that("directional LPs and vmax are partitioned by target block", {
  feasibility <- paste(
    deparse(body(RegCompassR:::.rc_directional_feasibility)),
    collapse = "\n"
  )
  vmax <- paste(
    deparse(body(RegCompassR:::.rc_build_microcompass_vmax_cache)),
    collapse = "\n"
  )
  expect_match(feasibility, ".rc_lp_block_size", fixed = TRUE)
  expect_match(feasibility, "rc_parallel_lapply", fixed = TRUE)
  expect_match(vmax, "selected_rows", fixed = TRUE)
  expect_match(vmax, "row_ids", fixed = TRUE)
  expect_match(vmax, "blocks = length(tasks)", fixed = TRUE)
})

test_that("Layer 2 propagates its monitor into engine and controls", {
  body_text <- paste(deparse(body(rc_regcompass_step_layer2)), collapse = "\n")
  expect_match(body_text, ".rc_run_microcompass_monitored", fixed = TRUE)
  expect_match(body_text, ".rc_run_microcompass_engine_monitored", fixed = TRUE)
  expect_match(body_text, ".rc_compact_meta_modules_for_layer2", fixed = TRUE)
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
