#' Run RegCompass from species-aware defaults
#'
#' Loads a Human-GEM or Mouse-GEM when `gem` is omitted, prepares the requested
#' medium, and calls [rc_run_regcompass()]. Pando routing is resolved separately
#' for every retained broad cell type. Stage 2 builds one multimodal WNN graph
#' per broad cell type and keeps final metacells condition-pure. Layer 2 uses
#' original MATLAB CORDA2 by default; FASTCORE and full-GEM scoring remain
#' explicit supplementary routes through arguments passed in `...`.
#'
#' When both medium arguments are omitted, Human-GEM uses
#' `"normal_human_plasma"` and Mouse-GEM uses `"mouse_plasma"`. Users may pass
#' any supported biological scenario or a prebuilt custom medium table.
#'
#' @param object A paired-cell Seurat RNA+ATAC object.
#' @param outdir Persistent output directory.
#' @param genome Genome object matching the selected species, ATAC coordinates,
#'   and regulatory regions.
#' @param species `"human"` or `"mouse"`.
#' @param gem Optional prepared species GEM.
#' @param gem_version Model release used when `gem` is omitted.
#' @param gem_source GEM source: automatic, bundled-only, or download.
#' @param pfm Optional motif collection. Pando's bundled motifs are used when
#'   omitted.
#' @param fragment_files Must be `FALSE`; metacell RNA and ATAC counts are
#'   aggregated from the existing assays.
#' @param medium_scenario Optional built-in biological scenario.
#' @param medium_scenarios Optional prebuilt user-defined or built-in medium table.
#' @param progress Show stage and total progress.
#' @param ... Arguments passed to [rc_run_regcompass()], including explicit
#'   supplementary Layer 2 mode controls when required.
#' @return A complete RegCompass result.
#' @export
rc_run_regcompass_one_shot <- function(
    object, outdir, genome,
    species = c("human", "mouse"),
    gem = NULL,
    gem_version = NULL,
    gem_source = c("auto", "bundled", "download"),
    pfm = NULL,
    fragment_files = FALSE,
    medium_scenario = NULL,
    medium_scenarios = NULL,
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
    species <- .rc_infer_gem_species(gem, species)
  }
  if (is.null(medium_scenarios)) {
    if (is.null(medium_scenario)) {
      medium_scenario <- if (identical(species, "human")) {
        "normal_human_plasma"
      } else {
        "mouse_plasma"
      }
    }
    medium_scenarios <- rc_make_medium_scenarios(
      gem = gem,
      scenario = medium_scenario,
      species = species
    )
  }
  rc_run_regcompass(
    object = object,
    gem = gem,
    outdir = outdir,
    genome = genome,
    pfm = pfm,
    species = species,
    fragment_files = fragment_files,
    medium_scenarios = medium_scenarios,
    progress = progress,
    ...
  )
}
