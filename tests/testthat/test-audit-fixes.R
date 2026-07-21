test_that("missing and explicit zero expression have identical strict penalties", {
  expression <- matrix(
    c(NA_real_, 0, 0.5, 4),
    nrow = 4,
    dimnames = list(c("missing", "zero", "low", "high"), "u1")
  )
  out <- rc_compute_multiome_penalty(expression)

  expect_equal(out$penalty["missing", "u1"], 1)
  expect_equal(out$penalty["zero", "u1"], 1)
  expect_equal(
    out$components$effective_reaction_expression["missing", "u1"],
    0
  )
  expect_true(out$components$missing_expression_flag["missing", "u1"])
  expect_false(out$components$missing_expression_flag["zero", "u1"])
  expect_true(all(out$penalty[c("low", "high"), "u1"] < 1))
  expect_error(
    rc_compute_multiome_penalty(expression, missing_penalty = 2),
    "must remain 1"
  )
})

test_that("GPR diagnostics require at least one complete isozyme group", {
  gpr <- list(R1 = list(c("g1", "g2"), c("g3", "g4")))

  incomplete <- rc_gpr_diagnostics(gpr, c("g1", "g3"))
  expect_equal(incomplete$n_complete_and_groups, 0)
  expect_true(incomplete$capacity_missing_flag)
  expect_true(incomplete$incomplete_and_group_flag)

  complete <- rc_gpr_diagnostics(gpr, c("g1", "g2"))
  expect_equal(complete$n_complete_and_groups, 1)
  expect_false(complete$capacity_missing_flag)
  expect_true(complete$incomplete_and_group_flag)
})

test_that("parallel helper rejects logical TRUE as BPPARAM", {
  expect_error(
    rc_parallel_lapply(1:2, identity, BPPARAM = TRUE),
    "logical TRUE is not valid"
  )
  expect_equal(
    unlist(rc_parallel_lapply(1:2, identity, BPPARAM = FALSE)),
    1:2
  )
})

test_that("zero-count ATAC features are removed by default", {
  rna <- Matrix::Matrix(
    matrix(c(1, 2, 3, 4), nrow = 2,
           dimnames = list(c("g1", "g2"), c("c1", "c2"))),
    sparse = TRUE
  )
  atac <- Matrix::Matrix(
    matrix(c(0, 0, 1, 2), nrow = 2, byrow = TRUE,
           dimnames = list(c("p0", "p1"), c("c1", "c2"))),
    sparse = TRUE
  )
  object <- SeuratObject::CreateSeuratObject(counts = rna)
  object[["ATAC"]] <- SeuratObject::CreateAssayObject(counts = atac)

  filtered <- .rc_drop_zero_count_atac_features(object, "ATAC")
  expect_identical(rownames(filtered$object[["ATAC"]]), "p1")
  expect_equal(filtered$diagnostics$n_input_peaks, 2)
  expect_equal(filtered$diagnostics$n_zero_count_peaks_excluded, 1)
  expect_equal(filtered$diagnostics$n_retained_peaks, 1)
})

test_that("condition-only pool identifiers do not depend on sample labels", {
  cells <- paste0("c", seq_len(6))
  counts <- Matrix::Matrix(
    matrix(seq_len(12), nrow = 2,
           dimnames = list(c("g1", "g2"), cells)),
    sparse = TRUE
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object$sample_id <- c("s1", "s1", "s2", "s3", "s4", "s4")
  object$condition <- c("A", "A", "A", "B", "B", "B")
  object$cell_type <- "T"

  first <- .rc_prepare_condition_only_object(object, "condition")
  object$sample_id <- rev(object$sample_id)
  second <- .rc_prepare_condition_only_object(object, "condition")

  expect_identical(first$sample_col, second$sample_col)
  expect_identical(
    first$object@meta.data[[first$sample_col]],
    second$object@meta.data[[second$sample_col]]
  )
  expect_setequal(
    unique(first$object@meta.data[[first$sample_col]]),
    c("A__condition_pool", "B__condition_pool")
  )
})

test_that("canonical metacell APIs make sample metadata optional", {
  expect_null(eval(formals(rc_regcompass_step_metacells)$sample_col))
  expect_null(eval(formals(rc_run_regcompass)$sample_col))

  wrapper_text <- paste(
    deparse(body(.rc_make_condition_pooled_metacells)),
    collapse = "\n"
  )
  expect_match(wrapper_text, ".rc_prepare_condition_only_object", fixed = TRUE)
  expect_match(wrapper_text, "ignored_balance_args", fixed = TRUE)
  expect_false(grepl(
    ".rc_balance_condition_celltype_cells",
    wrapper_text,
    fixed = TRUE
  ))
})

test_that("hidden inference-unit option is retired", {
  body_text <- paste(deparse(body(rc_run_microcompass)), collapse = "\n")
  expect_match(body_text, "retired `RegCompassR.inference_unit` option is ignored",
               fixed = TRUE)
  expect_match(body_text, "options(RegCompassR.inference_unit = NULL)",
               fixed = TRUE)
})

test_that("stepwise stages expose the requested parallel contract", {
  for (fun in list(
    rc_regcompass_step_meta_modules,
    rc_regcompass_step_layer1,
    rc_regcompass_step_layer2
  )) {
    f <- formals(fun)
    expect_identical(eval(f$parallel), TRUE)
    expect_null(eval(f$BPPARAM))
  }
})
