test_that("supplied Pando resources are preserved before dispatch", {
  pfm <- structure(list(name = "pfm"), class = "test_pfm")
  regions <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1L, end = 10L)
  )
  motif_tfs <- data.frame(
    motif = "M1", tf = "TF1", stringsAsFactors = FALSE
  )

  resources <- .rc_materialize_pando_resources(
    pfm = pfm,
    species = "human",
    pando_initiate_args = list(
      exclude_exons = TRUE,
      regions = regions
    ),
    pando_motif_args = list(
      motif_tfs = motif_tfs,
      cache_dir = "cache"
    )
  )

  expect_identical(resources$pfm, pfm)
  expect_identical(resources$pando_initiate_args$regions, regions)
  expect_identical(resources$pando_initiate_args$exclude_exons, TRUE)
  expect_identical(resources$pando_motif_args$motif_tfs, motif_tfs)
  expect_identical(resources$pando_motif_args$cache_dir, "cache")
})

test_that("Pando data loading never dereferences namespace lazydata", {
  loader_source <- paste(deparse(body(.rc_pando_data_object)), collapse = "\n")
  expect_false(grepl("lazydata", tolower(loader_source), fixed = TRUE))
  expect_true(grepl("utils::data", loader_source, fixed = TRUE))
  expect_true(grepl("lib.loc", loader_source, fixed = TRUE))
})

test_that("human Pando default regions are the phastCons plus SCREEN union", {
  region_source <- paste(
    deparse(body(.rc_default_pando_regions)), collapse = "\n"
  )
  expect_match(
    region_source,
    "phastConsElements20Mammals.UCSC.hg38",
    fixed = TRUE
  )
  expect_match(region_source, "SCREEN.ccRE.UCSC.hg38", fixed = TRUE)
  expect_match(region_source, "BiocGenerics::union", fixed = TRUE)
})

test_that("Stage 1 materializes resources before Pando job dispatch", {
  stage_source <- paste(
    deparse(body(rc_regcompass_step_grn)), collapse = "\n"
  )
  materialize <- regexpr(
    ".rc_materialize_pando_resources", stage_source, fixed = TRUE
  )[[1L]]
  dispatch <- regexpr(
    ".rc_fit_pando_by_celltype_route", stage_source, fixed = TRUE
  )[[1L]]
  expect_gt(materialize, 0L)
  expect_gt(dispatch, materialize)
  expect_true(grepl(
    "dispatch_extra_args$pando_initiate_args <- resources$pando_initiate_args",
    stage_source,
    fixed = TRUE
  ))
  expect_true(grepl(
    "extra_args = dispatch_extra_args",
    stage_source,
    fixed = TRUE
  ))
  expect_true(grepl(
    "stored_in_step_params = FALSE",
    stage_source,
    fixed = TRUE
  ))
})

test_that("standard and condition Pando routes receive the same materialized regions", {
  route_source <- paste(
    deparse(body(.rc_fit_pando_by_celltype_route)), collapse = "\n"
  )
  expect_match(route_source, ".rc_run_condition_pando_batch", fixed = TRUE)
  expect_match(route_source, ".rc_run_standard_pando_celltype_job", fixed = TRUE)
  occurrences <- gregexpr(
    "extra_args = extra_args", route_source, fixed = TRUE
  )[[1L]]
  occurrences <- occurrences[occurrences > 0L]
  expect_gte(length(occurrences), 2L)

  condition_source <- paste(
    deparse(body(.rc_run_condition_pando_batch)), collapse = "\n"
  )
  standard_source <- paste(
    deparse(body(.rc_run_standard_pando_celltype_job)), collapse = "\n"
  )
  expect_match(
    condition_source,
    "args <- c(base[setdiff(names(base), names(extra_args))], extra_args)",
    fixed = TRUE
  )
  expect_match(
    standard_source,
    "args <- c(base[setdiff(names(base), names(job_extra))], job_extra)",
    fixed = TRUE
  )
})
