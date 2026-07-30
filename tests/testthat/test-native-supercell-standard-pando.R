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
  design <- RegCompassR:::.rc_resolve_condition_design(object, "condition")
  expect_identical(design$analysis_mode, "standard_pando")
  expect_identical(design$fallback_reason, "condition_col_absent")
  expect_true(design$condition_col %in% colnames(design$object@meta.data))
})

test_that("native SuperCell receives separate condition and cell-type inputs", {
  fun <- get(
    ".rc_native_supercell_membership",
    envir = asNamespace("RegCompassR"),
    inherits = FALSE
  )
  text <- paste(deparse(body(fun)), collapse = "\n")
  expect_match(text, "SCimplify_from_embedding", fixed = TRUE)
  expect_match(text, "cell.annotation", fixed = TRUE)
  expect_match(text, "cell.split.condition", fixed = TRUE)
  expect_false(grepl("condition__cell_type", text, fixed = TRUE))
  expect_false(grepl("stratum_col", text, fixed = TRUE))
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
