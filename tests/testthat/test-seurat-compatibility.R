test_that("Seurat stack profiles retain v4 default and accept v5", {
  v4 <- .rc_validate_seurat_stack_versions(c(
    SeuratObject = "4.1.4",
    Seurat = "4.4.0",
    Signac = "1.11.0"
  ))
  expect_identical(v4$profile, "seurat_v4_default")
  expect_true(v4$default)

  v4_later <- .rc_validate_seurat_stack_versions(c(
    SeuratObject = "4.1.5",
    Seurat = "4.4.1",
    Signac = "1.13.0"
  ))
  expect_identical(v4_later$profile, "seurat_v4_default")

  v5 <- .rc_validate_seurat_stack_versions(c(
    SeuratObject = "5.0.2",
    Seurat = "5.0.3",
    Signac = "1.14.0"
  ))
  expect_identical(v5$profile, "seurat_v5_compatible")
  expect_false(v5$default)

  expect_error(
    .rc_validate_seurat_stack_versions(c(
      SeuratObject = "5.0.2",
      Seurat = "4.4.0",
      Signac = "1.12.0"
    )),
    "same major version"
  )
  expect_error(
    .rc_validate_seurat_stack_versions(c(
      SeuratObject = "5.0.2",
      Seurat = "5.0.3",
      Signac = "1.11.0"
    )),
    "Signac >=1.12.0"
  )
  expect_error(
    .rc_validate_seurat_stack_versions(c(
      SeuratObject = "5.0.2",
      Seurat = "5.0.3",
      Signac = "2.0.0"
    )),
    "ChromatinAssay5"
  )
  unsupported <- .rc_validate_seurat_stack_versions(c(
    SeuratObject = "6.0.0",
    Seurat = "6.0.0",
    Signac = "1.17.0"
  ), error = FALSE)
  expect_identical(unsupported$profile, "unsupported")
  expect_false(unsupported$supported)
})

test_that("v3 Assay access is stable under either Seurat major", {
  old_option <- getOption("Seurat.object.assay.version")
  on.exit(options(Seurat.object.assay.version = old_option), add = TRUE)
  options(Seurat.object.assay.version = "v3")

  counts <- Matrix::Matrix(
    matrix(c(1, 0, 2, 3, 0, 1, 4, 0, 2, 1, 0, 3), nrow = 3),
    sparse = TRUE,
    dimnames = list(
      paste0("gene", 1:3),
      paste0("cell", 1:4)
    )
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  expect_false(inherits(object[["RNA"]], "Assay5"))
  expect_equal(
    as.matrix(.rc_get_assay_counts(object, "RNA")),
    as.matrix(counts)
  )

  normalized <- log1p(counts)
  object <- .rc_set_assay_matrix(
    object, assay = "RNA", layer = "data", new_data = normalized
  )
  expect_equal(
    as.matrix(.rc_get_assay_matrix(object, "RNA", "data")),
    as.matrix(normalized)
  )

  object <- .rc_prepare_seurat_assays(
    object,
    assays = "RNA",
    required_layers = "counts",
    optional_layers = "data"
  )
  provenance <- object@misc$regcompass_seurat_compatibility
  expect_identical(provenance$default_input_profile, "Seurat_v4_v3_Assay")
  expect_identical(unname(provenance$assay_storage[["RNA"]]), "v3_slots")
})

test_that("Assay5 split layers are joined on a working object", {
  create_assay5 <- .rc_seurat_object_export("CreateAssay5Object")
  layers_fun <- .rc_seurat_object_export("Layers")
  skip_if(is.null(create_assay5) || is.null(layers_fun),
          "SeuratObject v5 Assay5 API is unavailable")

  old_option <- getOption("Seurat.object.assay.version")
  on.exit(options(Seurat.object.assay.version = old_option), add = TRUE)
  options(Seurat.object.assay.version = "v5")

  counts <- Matrix::Matrix(
    matrix(c(1, 0, 2, 3, 0, 1, 4, 0, 2, 1, 0, 3), nrow = 3),
    sparse = TRUE,
    dimnames = list(
      paste0("gene", 1:3),
      paste0("cell", 1:4)
    )
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  expect_true(inherits(object[["RNA"]], "Assay5"))
  batch <- factor(c("A", "A", "B", "B"))
  object[["RNA"]] <- split(object[["RNA"]], f = batch)
  before <- .rc_assay_layer_names(object[["RNA"]])
  expect_gt(sum(startsWith(before, "counts.")), 1L)

  prepared <- .rc_prepare_seurat_assays(
    object,
    assays = "RNA",
    required_layers = "counts",
    optional_layers = "data"
  )
  after <- .rc_assay_layer_names(prepared[["RNA"]])
  expect_identical(after[after == "counts"], "counts")
  expect_false(any(startsWith(after, "counts.")))
  expect_equal(
    as.matrix(.rc_get_assay_counts(prepared, "RNA")),
    as.matrix(counts[, colnames(prepared), drop = FALSE])
  )

  normalized <- log1p(.rc_get_assay_counts(prepared, "RNA"))
  prepared <- .rc_set_assay_matrix(
    prepared, assay = "RNA", layer = "data", new_data = normalized
  )
  expect_equal(
    as.matrix(.rc_get_assay_matrix(prepared, "RNA", "data")),
    as.matrix(normalized)
  )
  provenance <- prepared@misc$regcompass_seurat_compatibility
  expect_identical(unname(provenance$assay_storage[["RNA"]]), "v5_layers")
  expect_true("RNA:counts" %in% names(provenance$joined_layers))

  original_layers <- .rc_assay_layer_names(object[["RNA"]])
  expect_gt(sum(startsWith(original_layers, "counts.")), 1L)
})

test_that("DESCRIPTION and README preserve the default profile contract", {
  description <- utils::packageDescription("RegCompassR")
  expect_identical(
    description[["Config/RegCompassR/DefaultSeuratProfile"]],
    "seurat_v4_default"
  )
  expect_identical(
    description[["Config/RegCompassR/SupportedSeuratMajors"]],
    "4,5"
  )
  expect_identical(
    description[["Config/RegCompassR/SeuratV5MinSignacVersion"]],
    "1.12.0"
  )

  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  roots <- unique(c(
    if (nzchar(workspace)) workspace else character(),
    ".", "..", file.path("..", "..")
  ))
  roots <- roots[file.exists(file.path(roots, "README.md"))]
  if (!length(roots)) skip("Source README is unavailable")
  readme <- paste(readLines(file.path(roots[[1L]], "README.md"), warn = FALSE),
                  collapse = "\n")
  expect_match(readme, "Default validated profile: Seurat v4", fixed = TRUE)
  expect_match(readme, "Optional compatible profile: Seurat v5", fixed = TRUE)
  expect_match(readme, 'Seurat.object.assay.version = "v3"', fixed = TRUE)
  expect_match(readme, "Signac 2.x `ChromatinAssay5` is not yet supported",
               fixed = TRUE)
})
