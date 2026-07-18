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

test_that("species arguments preserve legacy positional ordering", {
  workflow_formals <- names(formals(rc_run_regcompass))
  one_shot_formals <- names(formals(rc_run_regcompass_one_shot))

  expect_lt(
    match("sample_col", workflow_formals),
    match("species", workflow_formals)
  )
  expect_lt(match("gem", one_shot_formals), match("species", one_shot_formals))
  expect_lt(
    match("medium_scenarios", one_shot_formals),
    match("species", one_shot_formals)
  )
})

test_that("stratum confidence forwards the selected RNA assay", {
  text <- paste(deparse(body(.rc_run_regcompass_stratum)), collapse = "\n")
  expect_match(
    text,
    "atac_assay = atac_assay,\\s*rna_assay = rna_assay"
  )
})

test_that("multi-sample projected edges receive sample-local module IDs", {
  input <- data.frame(
    sample_id = c("s1", "s1", "s2", "s2"),
    tf = c("TF1", "TF1", "TF2", "TF2"),
    target = c("A", "B", "A", "C"),
    estimate = c(1, 1, 1, 1),
    stringsAsFactors = FALSE
  )
  projected <- rc_project_metabolic_grn(
    input,
    metabolic_genes = c("A", "B", "C"),
    top_k = 5,
    min_shared_tfs = 1,
    min_tf_jaccard = 0,
    include_direct_metabolic_tf = FALSE
  )
  expect_true(all(
    startsWith(
      projected$edges$module_id,
      paste0(projected$edges$sample_id, "::")
    )
  ))
})

test_that("signed eligibility is applied before top-k component pruning", {
  input <- data.frame(
    sample_id = rep("s1", 6),
    tf = c("TF1", "TF1", "TF2", "TF2", "TF3", "TF3"),
    target = c("A", "B", "A", "C", "C", "D"),
    estimate = c(10, -10, 1, 1, 10, -10),
    stringsAsFactors = FALSE
  )
  projected <- rc_project_metabolic_grn(
    input,
    metabolic_genes = c("A", "B", "C", "D"),
    top_k = 1,
    min_shared_tfs = 1,
    min_tf_jaccard = 0,
    include_direct_metabolic_tf = FALSE
  )
  pair <- paste(projected$edges$gene_a, projected$edges$gene_b, sep = "-")
  concordant <- pair == "A-C"
  discordant <- pair %in% c("A-B", "C-D")

  expect_true(any(concordant))
  expect_true(projected$edges$used_for_component[concordant])
  expect_true(all(!projected$edges$used_for_component[discordant]))
  expect_equal(
    projected$nodes$module_id[projected$nodes$gene == "A"],
    projected$nodes$module_id[projected$nodes$gene == "C"]
  )
})
