test_that("condition-metacell checkpoints require an identical cache contract", {
  outdir <- tempfile("regcompass-metacell-cache-")
  dir.create(file.path(outdir, "condition=A"), recursive = TRUE)
  file.create(file.path(
    outdir, "condition=A", "metacell_metadata.tsv.gz"
  ))
  file.create(file.path(outdir, "condition=A", "rna_counts.rds"))
  file.create(file.path(outdir, "condition=A", "atac_counts.rds"))

  contract <- list(
    schema_version = "regcompass_condition_metacell_cache_v1",
    condition_col = "condition",
    celltype_col = "cell_type",
    analysis_args = list(gamma = 30L)
  )

  expect_error(
    .rc_validate_condition_metacell_cache(
      outdir, contract, overwrite = FALSE
    ),
    "predate the audited cache contract"
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
    "different cells, condition/cell-type labels, assay contents"
  )
  expect_false(.rc_validate_condition_metacell_cache(
    outdir, changed, overwrite = TRUE
  ))
})

test_that("downstream stages reject legacy metacell provenance", {
  params <- list(
    condition_col = "condition",
    celltype_col = "cell_type",
    rna_assay = "RNA",
    atac_assay = "ATAC",
    metacell_args = list(gamma = 30L)
  )
  current <- structure(
    list(
      pooled = list(
        input_design = list(
          condition_only_stratification = TRUE,
          supercell_label_col = "cell_type",
          celltype_assignment = paste(
            "SuperCell2 label-guided construction followed by",
            "dominant membership assignment"
          ),
          gamma = 30L
        ),
        cache_contract = list(
          schema_version = "regcompass_condition_metacell_cache_v1",
          condition_col = "condition",
          celltype_col = "cell_type",
          rna_assay = "RNA",
          atac_assay = "ATAC",
          label_col = "cell_type",
          analysis_args = list(gamma = 30L)
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
    "legacy or incompatible metacell artifact"
  )
})

test_that("DESCRIPTION exact pins match the runtime version gate", {
  imports <- packageDescription("RegCompassR")$Imports
  expect_match(imports, "SeuratObject \\(= 4.1.4\\)")
  expect_match(imports, "Seurat \\(= 4.4.0\\)")
  expect_match(imports, "Signac \\(= 1.11.0\\)")
})

test_that("the public function index identifies the true stepwise tutorial", {
  path <- testthat::test_path("..", "..", "docs", "functions.md")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(text, "Level 2: true stepwise workflow", fixed = TRUE)
  expect_false(grepl("Level 2: saved-stage audit", text, fixed = TRUE))
})
