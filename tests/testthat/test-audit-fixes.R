test_that("missing and explicit zero expression have identical strict penalties", {
  expression <- matrix(c(NA_real_, 0, 0.5, 4), nrow = 4,
    dimnames = list(c("missing", "zero", "low", "high"), "u1"))
  out <- rc_compute_multiome_penalty(expression)
  expect_equal(out$penalty["missing", "u1"], 1)
  expect_equal(out$penalty["zero", "u1"], 1)
  expect_equal(out$components$effective_reaction_expression["missing", "u1"], 0)
  expect_true(out$components$missing_expression_flag["missing", "u1"])
  expect_false(out$components$missing_expression_flag["zero", "u1"])
  expect_true(all(out$penalty[c("low", "high"), "u1"] < 1))
  expect_error(rc_compute_multiome_penalty(expression, missing_penalty = 2), "must remain 1")
})

test_that("GPR diagnostics require at least one complete isozyme group", {
  gpr <- list(R1 = list(c("g1", "g2"), c("g3", "g4")))
  incomplete <- rc_gpr_diagnostics(gpr, c("g1", "g3"))
  expect_equal(incomplete$n_complete_and_groups, 0)
  expect_true(incomplete$capacity_missing_flag)
  complete <- rc_gpr_diagnostics(gpr, c("g1", "g2"))
  expect_equal(complete$n_complete_and_groups, 1)
  expect_false(complete$capacity_missing_flag)
})

test_that("parallel helper rejects logical TRUE as BPPARAM", {
  expect_error(rc_parallel_lapply(1:2, identity, BPPARAM = TRUE), "logical TRUE is not valid")
  expect_equal(unlist(rc_parallel_lapply(1:2, identity, BPPARAM = FALSE)), 1:2)
})

test_that("zero-count ATAC features are removed by default", {
  rna <- Matrix::Matrix(matrix(c(1, 2, 3, 4), nrow = 2,
    dimnames = list(c("g1", "g2"), c("c1", "c2"))), sparse = TRUE)
  atac <- Matrix::Matrix(matrix(c(0, 0, 1, 2), nrow = 2, byrow = TRUE,
    dimnames = list(c("p0", "p1"), c("c1", "c2"))), sparse = TRUE)
  object <- SeuratObject::CreateSeuratObject(counts = rna)
  object[["ATAC"]] <- SeuratObject::CreateAssayObject(counts = atac)
  filtered <- .rc_drop_zero_count_atac_features(object, "ATAC")
  expect_identical(rownames(filtered$object[["ATAC"]]), "p1")
  expect_equal(filtered$diagnostics$n_input_peaks, 2)
  expect_equal(filtered$diagnostics$n_zero_count_peaks_excluded, 1)
  expect_equal(filtered$diagnostics$n_retained_peaks, 1)
})

test_that("condition pool identifiers do not use sample labels", {
  first <- data.frame(sample_id = c("s1", "s2"), condition = c("A", "B"))
  second <- first
  second$sample_id <- rev(second$sample_id)
  expect_identical(.rc_condition_only_sample_col(first), .rc_condition_only_sample_col(second))
})

test_that("canonical metacell APIs exclude sample balancing", {
  expect_null(eval(formals(rc_regcompass_step_metacells)$sample_col))
  expect_null(eval(formals(rc_run_regcompass)$sample_col))
  wrapper_text <- paste(deparse(body(.rc_make_condition_pooled_metacells)), collapse = "\n")
  expect_match(wrapper_text, ".rc_condition_only_sample_col", fixed = TRUE)
  expect_match(wrapper_text, "Sample balancing is not part", fixed = TRUE)
  expect_false(grepl(".rc_balance_condition_celltype_cells", wrapper_text, fixed = TRUE))
})

test_that("hidden inference-unit option is retired", {
  body_text <- paste(deparse(body(rc_run_microcompass)), collapse = "\n")
  expect_match(body_text, "retired `RegCompassR.inference_unit` option is ignored", fixed = TRUE)
  expect_match(body_text, "options(RegCompassR.inference_unit = NULL)", fixed = TRUE)
})

test_that("parallel contracts apply to GRN, Layer 1 and Layer 2", {
  for (fun in list(rc_regcompass_step_grn, rc_regcompass_step_layer1, rc_regcompass_step_layer2)) {
    f <- formals(fun)
    expect_identical(eval(f$parallel), TRUE)
    expect_null(eval(f$BPPARAM))
  }
  expect_false(any(c("parallel", "BPPARAM") %in% names(formals(rc_regcompass_step_meta_modules))))
})
