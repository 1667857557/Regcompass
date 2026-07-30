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

test_that("native SuperCell builds independent cell-type graphs jointly across conditions", {
  fun <- get(
    ".rc_native_supercell_membership",
    envir = asNamespace("RegCompassR"),
    inherits = FALSE
  )
  text <- paste(deparse(body(fun)), collapse = "\n")
  expect_match(text, "SCimplify_by_graph_group_from_embedding", fixed = TRUE)
  expect_match(text, "cell.graph.group", fixed = TRUE)
  expect_match(text, "cell.split.condition", fixed = TRUE)
  expect_match(text, ".rc_scale_embedding_block_by_group", fixed = TRUE)
  expect_false(grepl("cell.annotation", text, fixed = TRUE))
  expect_false(grepl("condition__cell_type", text, fixed = TRUE))
  expect_false(grepl("stratum_col", text, fixed = TRUE))
})

test_that("cell-type scaling pools conditions but isolates other cell types", {
  x <- matrix(
    c(
      0, 1,
      2, 3,
      4, 5,
      100, 101,
      102, 103,
      104, 105
    ),
    ncol = 2,
    byrow = TRUE
  )
  group <- c("A", "A", "A", "B", "B", "B")
  first <- .rc_scale_embedding_block_by_group(x, group)
  changed <- x
  changed[group == "B", ] <- changed[group == "B", ] * 1000 + 5000
  second <- .rc_scale_embedding_block_by_group(changed, group)
  expect_equal(first[group == "A", ], second[group == "A", ], tolerance = 1e-12)
  expect_equal(colMeans(first[group == "A", , drop = FALSE]), c(0, 0), tolerance = 1e-12)
  expect_equal(colMeans(first[group == "B", , drop = FALSE]), c(0, 0), tolerance = 1e-12)
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
