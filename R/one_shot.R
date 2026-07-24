#' Run RegCompass from species-aware defaults
#'
#' Uses the bundled pinned Human-GEM or Mouse-GEM by default when `gem` is
#' omitted. Set `gem_source = "download"` when rebuilding from an upstream
#' release.
#'
#' @param object A Seurat RNA+ATAC object.
#' @param outdir Persistent output directory.
#' @param pfm Motif position-frequency matrices.
#' @param genome Genome object matching the selected species and ATAC coordinates.
#' @param fragment_files Must be `FALSE` for the canonical peak-count path.
#' @param species `"human"` or `"mouse"`.
#' @param gem Optional prebuilt species GEM.
#' @param gem_version Pinned model release.
#' @param gem_source GEM source: automatic, bundled-only, or download.
#' @param medium_scenario Medium preset identifier.
#' @param medium_scenarios Optional prebuilt medium table.
#' @param progress Show stage and total progress.
#' @param ... Arguments passed to [rc_run_regcompass()].
#' @return A canonical RegCompass result list.
#' @export
rc_run_regcompass_one_shot <- function(
    object, outdir, pfm, genome,
    fragment_files = FALSE,
    gem = NULL,
    medium_scenario = "physiologic",
    medium_scenarios = NULL,
    species = c("human", "mouse"),
    gem_version = NULL,
    gem_source = c("auto", "bundled", "download"),
    progress = getOption("RegCompassR.progress", TRUE),
    ...) {
  species <- match.arg(species)
  gem_source <- match.arg(gem_source)
  if (is.null(gem_version)) {
    gem_version <- if (identical(species, "human")) "2.0.0" else "1.8.0"
  }
  if (is.null(gem)) {
    gem <- rc_prepare_gem(
      species = species,
      version = gem_version,
      source = gem_source
    )
  } else {
    .rc_infer_gem_species(gem, species)
  }
  if (is.null(medium_scenarios)) {
    medium_scenarios <- rc_make_medium_scenarios(
      gem = gem, scenario = medium_scenario, species = species
    )
  }
  rc_run_regcompass(
    object = object, gem = gem, outdir = outdir,
    pfm = pfm, genome = genome,
    fragment_files = fragment_files,
    species = species,
    medium_scenarios = medium_scenarios,
    progress = progress,
    ...
  )
}
