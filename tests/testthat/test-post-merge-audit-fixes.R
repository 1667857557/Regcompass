test_that("full-GEM cache identity changes with GEM structure", {
  make_gem <- function(stoichiometry, version) {
    S <- Matrix::Matrix(
      matrix(c(-1, stoichiometry), nrow = 1),
      sparse = TRUE,
      dimnames = list("m_e", c("EX_m", "R1"))
    )
    list(
      S = S,
      lb = stats::setNames(c(-1000, 0), colnames(S)),
      ub = stats::setNames(c(1000, 1000), colnames(S)),
      reaction_meta = data.frame(
        reaction_id = colnames(S),
        role = c("exchange", "internal"),
        stringsAsFactors = FALSE
      ),
      model_info = list(
        species = "human",
        source = "test/GEM",
        version = version,
        commit = version,
        checksum = paste0("checksum-", version)
      )
    )
  }
  medium <- data.frame(
    medium_scenario_id = "base",
    exchange_reaction_id = NA_character_,
    lb = NA_real_,
    ub = NA_real_,
    available = FALSE,
    .no_constraints = TRUE,
    stringsAsFactors = FALSE
  )
  dirs <- data.frame(
    reaction_id = "R1",
    target_direction = "forward",
    stringsAsFactors = FALSE
  )
  cache_dir <- tempfile("full-gem-cache-")
  first <- rc_build_full_gem_cache(
    make_gem(1, "v1"), dirs, medium, cache_dir = cache_dir
  )
  second <- rc_build_full_gem_cache(
    make_gem(2, "v2"), dirs, medium, cache_dir = cache_dir
  )
  first_file <- attr(first, "summary")$file[[1L]]
  second_file <- attr(second, "summary")$file[[1L]]

  expect_false(identical(first_file, second_file))
  expect_false(identical(
    attr(first, "summary")$gem_fingerprint,
    attr(second, "summary")$gem_fingerprint
  ))
  expect_equal(as.numeric(readRDS(second_file)$S["m_e", "R1"]), 2)
})

test_that("incompatible legacy species caches are invalidated", {
  S <- Matrix::Matrix(
    matrix(c(-1, 1), nrow = 1),
    sparse = TRUE,
    dimnames = list("m_e", c("EX_m", "R1"))
  )
  legacy <- list(
    S = S,
    lb = stats::setNames(c(-1000, 0), colnames(S)),
    ub = stats::setNames(c(1000, 1000), colnames(S)),
    model_info = list(source = "SysBioChalmers/Human-GEM", version = "2.0.0")
  )
  file <- tempfile(fileext = ".rds")
  saveRDS(legacy, file)
  spec <- .rc_species_gem_spec("human", "2.0.0")

  expect_warning(
    cached <- .rc_load_compatible_species_gem(file, spec),
    "Removing incompatible cached"
  )
  expect_null(cached)
  expect_false(file.exists(file))
})

test_that("public runner arguments follow processing order", {
  expect_identical(
    names(formals(rc_run_regcompass)),
    c(
      "object", "gem", "outdir", "genome", "pfm", "species",
      "condition_col", "celltype_col", "rna_assay", "atac_assay",
      "pando_args",
      "sample_col", "fragment_files", "metacell_args",
      "meta_module_args",
      "layer1_args",
      "medium_scenarios", "model_mode", "layer2_args",
      "upstream_workers", "layer2_workers", "progress"
    )
  )
  expect_identical(
    names(formals(rc_run_regcompass_one_shot)),
    c(
      "object", "outdir", "genome",
      "species", "gem", "gem_version", "gem_source",
      "pfm", "fragment_files",
      "medium_scenario", "medium_scenarios",
      "progress", "..."
    )
  )
})
