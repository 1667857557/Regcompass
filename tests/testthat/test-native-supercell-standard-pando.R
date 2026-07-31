test_that("condition design selects the intended Pando mode", {
  counts <- Matrix::Matrix(
    matrix(
      c(1, 0, 2, 0, 1, 1),
      nrow = 2,
      dimnames = list(c("g1", "g2"), paste0("cell", 1:3))
    ),
    sparse = TRUE
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object$cell_type <- c("T", "T", "T")
  object$condition <- c("A", "A", "B")
  design <- RegCompassR:::.rc_resolve_condition_design(object, "condition")
  expect_identical(design$analysis_mode, "condition_grn")
  expect_setequal(design$condition_levels, c("A", "B"))

  object$condition <- "A"
  design <- RegCompassR:::.rc_resolve_condition_design(object, "condition")
  expect_identical(design$analysis_mode, "standard_pando")
  expect_identical(design$fallback_reason, "fewer_than_two_condition_levels")

  object@meta.data$condition <- NULL
  expect_error(
    RegCompassR:::.rc_resolve_condition_design(object, "condition"),
    "Explicitly requested `condition_col` was not found"
  )
  design <- RegCompassR:::.rc_resolve_condition_design(object, NULL)
  expect_identical(design$analysis_mode, "standard_pando")
  expect_identical(design$fallback_reason, "condition_col_not_supplied")
  expect_true(design$condition_col %in% colnames(design$object@meta.data))
})

test_that("metacell workflow delegates graphing and aggregation to upstream SuperCell2", {
  membership_fun <- get(
    ".rc_native_supercell_membership",
    envir = asNamespace("RegCompassR"),
    inherits = FALSE
  )
  aggregation_fun <- get(
    ".rc_aggregate_native_metacell_counts",
    envir = asNamespace("RegCompassR"),
    inherits = FALSE
  )
  membership_text <- paste(deparse(body(membership_fun)), collapse = "\n")
  aggregation_text <- paste(deparse(body(aggregation_fun)), collapse = "\n")
  expect_match(membership_text, "SCimplify_for_Seurat", fixed = TRUE)
  expect_match(membership_text, "return.seurat = FALSE", fixed = TRUE)
  expect_match(membership_text, "paste(parent[cells], condition", fixed = TRUE)
  expect_match(aggregation_text, "membership = numeric_membership", fixed = TRUE)
  expect_match(aggregation_text, "return.seurat = TRUE", fixed = TRUE)
  expect_false(grepl("SCimplify_by_graph_group_from_embedding", membership_text,
                     fixed = TRUE))
  expect_false(grepl(".rc_scale_embedding_block_by_group", membership_text,
                     fixed = TRUE))
})

test_that("native SuperCell2 defaults retain WNN-oriented upstream controls", {
  defaults <- RegCompassR:::.rc_condition_metacell_defaults()
  expect_identical(defaults$gamma, 20L)
  expect_identical(defaults$k.knn, 30L)
  expect_true(defaults$kernel)
  expect_null(defaults$kith)
  expect_false(defaults$metacellNormalization)
  expect_false(defaults$avg.in.data)
})

test_that("one-condition standard Pando fallback sanitizes condition arguments", {
  args <- RegCompassR:::.rc_standard_pando_infer_args(list(
    candidate_screen = "motif_domain",
    condition_mix = 0.5,
    condition_weight = "equal",
    outer_nfolds = 5L,
    inner_nfolds = 5L,
    lambda_selection = "lambda.1se",
    scale = TRUE,
    family = "gaussian"
  ))
  expect_false(args$scale)
  expect_identical(args$interaction_term, ":")
  expect_identical(args$family, "gaussian")
  expect_false(any(c(
    "candidate_screen", "condition_mix", "condition_weight",
    "outer_nfolds", "inner_nfolds", "lambda_selection"
  ) %in% names(args)))
  adjustment <- attr(args, "standard_fallback_adjustments")
  expect_true(adjustment$scale_forced_false)
  expect_setequal(
    adjustment$dropped_condition_arguments,
    c(
      "candidate_screen", "condition_mix", "condition_weight",
      "outer_nfolds", "inner_nfolds", "lambda_selection"
    )
  )
})

test_that("standard Pando path calculates no condition coefficients", {
  fit <- get(
    ".rc_fit_standard_pando_by_cell_type",
    envir = asNamespace("RegCompassR"),
    inherits = FALSE
  )
  text <- paste(deparse(body(fit)), collapse = "\n")
  expect_match(text, "Pando::infer_grn", fixed = TRUE)
  expect_match(text, "condition_coefficients_calculated = FALSE", fixed = TRUE)
  expect_false(grepl("Pando::infer_condition_grn", text, fixed = TRUE))
})

test_that("package no longer collates runtime override files", {
  root <- normalizePath(file.path("..", ".."), mustWork = TRUE)
  description <- paste(
    readLines(file.path(root, "DESCRIPTION"), warn = FALSE),
    collapse = "\n"
  )
  expect_false(grepl("zzz00_absolute_pando_contract", description, fixed = TRUE))
  expect_false(grepl("zzz01_fixed_gamma_metacells", description, fixed = TRUE))
  expect_false(grepl("zzz02_layer1_policy", description, fixed = TRUE))
  expect_false(grepl("zzz03_compass_gpr_penalty", description, fixed = TRUE))
  expect_false(grepl("zzz04_canonical_pando_fit_schema", description, fixed = TRUE))
})
