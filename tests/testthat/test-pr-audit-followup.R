test_that("cell-type graph checkpoints require an identical cache contract", {
  outdir <- tempfile("regcompass-metacell-cache-")
  dir.create(outdir, recursive = TRUE)
  file.create(file.path(outdir, "metacell_metadata.tsv.gz"))
  file.create(file.path(outdir, "membership.tsv.gz"))
  file.create(file.path(outdir, "rna_counts.rds"))
  file.create(file.path(outdir, "atac_counts.rds"))
  file.create(file.path(outdir, "metacell_object.rds"))

  contract <- list(
    schema_version = "regcompass_celltype_graph_condition_joint_cache_v2",
    condition_col = "condition",
    celltype_col = "cell_type",
    graph_scope = "one_independent_graph_per_cell_type",
    condition_scope = "all_conditions_joint_within_cell_type_graph",
    analysis_args = list(gamma = 30L)
  )
  expect_error(
    .rc_validate_condition_metacell_cache(
      outdir, contract, overwrite = FALSE
    ),
    "different cell-type graph contract"
  )
  saveRDS(
    contract,
    file.path(outdir, "condition_metacell_cache_contract.rds")
  )
  expect_true(.rc_validate_condition_metacell_cache(
    outdir, contract, overwrite = FALSE
  ))
  changed <- contract
  changed$analysis_args$gamma <- 31L
  expect_error(
    .rc_validate_condition_metacell_cache(
      outdir, changed, overwrite = FALSE
    ),
    "different cell-type graph contract"
  )
  expect_false(.rc_validate_condition_metacell_cache(
    outdir, changed, overwrite = TRUE
  ))
})

test_that("matrix fingerprints detect value changes beyond marginal sums", {
  first <- matrix(
    c(1, 0, 0, 1), nrow = 2,
    dimnames = list(c("g1", "g2"), c("u1", "u2"))
  )
  second <- matrix(
    c(0, 1, 1, 0), nrow = 2,
    dimnames = dimnames(first)
  )
  first_fingerprint <- .rc_condition_metacell_matrix_fingerprint(first)
  second_fingerprint <- .rc_condition_metacell_matrix_fingerprint(second)
  expect_identical(
    first_fingerprint$row_sums_md5,
    second_fingerprint$row_sums_md5
  )
  expect_identical(
    first_fingerprint$col_sums_md5,
    second_fingerprint$col_sums_md5
  )
  expect_false(identical(
    first_fingerprint$values_md5,
    second_fingerprint$values_md5
  ))
})

test_that("downstream stages require cell-type-independent condition-joint provenance", {
  params <- list(
    condition_col = "condition",
    celltype_col = "cell_type",
    rna_assay = "RNA",
    atac_assay = "ATAC",
    analysis_mode = "condition_grn",
    metacell_args = list(gamma = 30L)
  )
  current <- structure(
    list(
      pooled = list(
        input_design = list(
          native_supercell_api =
            "SCimplify_by_graph_group_from_embedding",
          graph_group_argument = "cell.graph.group",
          condition_argument = "cell.split.condition",
          graph_scope = "one_independent_graph_per_cell_type",
          condition_scope = "all_conditions_joint_within_cell_type_graph",
          membership_split_timing = "after_joint_graph_clustering",
          embedding_scaling =
            "within_celltype_joint_condition_equal_modality_blocks",
          temporary_combined_stratum = FALSE
        ),
        cache_contract = list(
          schema_version =
            "regcompass_celltype_graph_condition_joint_cache_v2",
          condition_col = "condition",
          celltype_col = "cell_type",
          rna_assay = "RNA",
          atac_assay = "ATAC"
        )
      ),
      params = params
    ),
    class = c("regcompass_metacell_step", "list")
  )
  expect_invisible(.rc_require_stage_class(
    current,
    "regcompass_metacell_step",
    "metacells",
    "rc_regcompass_step_metacells"
  ))
  legacy <- current
  legacy$pooled$cache_contract <- NULL
  expect_error(
    .rc_require_stage_class(
      legacy,
      "regcompass_metacell_step",
      "metacells",
      "rc_regcompass_step_metacells"
    ),
    "not a cell-type-independent, condition-joint SuperCell artifact"
  )
})

test_that("DESCRIPTION preserves the default v4 profile and supported majors", {
  description <- packageDescription("RegCompassR")
  imports <- description$Imports
  expect_match(imports, "SeuratObject \\(>= 4.1.4\\)")
  expect_match(imports, "Seurat \\(>= 4.4.0\\)")
  expect_match(imports, "Signac \\(>= 1.11.0\\)")
  expect_identical(
    description[["Config/RegCompassR/DefaultSeuratProfile"]],
    "seurat_v4_default"
  )
  expect_identical(
    description[["Config/RegCompassR/SupportedSeuratMajors"]],
    "4,5"
  )
  expect_identical(
    description[["Config/RegCompassR/SupportedSignacMajor"]],
    "1"
  )
  expect_identical(
    description[["Config/RegCompassR/SeuratObjectVersion"]],
    "4.1.4"
  )
  expect_identical(
    description[["Config/RegCompassR/SeuratVersion"]],
    "4.4.0"
  )
  expect_identical(
    description[["Config/RegCompassR/SignacVersion"]],
    "1.11.0"
  )
})

test_that("the public function index identifies automatic and stepwise modes", {
  path <- testthat::test_path("..", "..", "docs", "functions.md")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(text, "Complete workflows", fixed = TRUE)
  expect_match(text, "Restartable stages", fixed = TRUE)
  expect_match(text, "standard_pando", fixed = TRUE)
  expect_match(text, "condition_grn", fixed = TRUE)
})
