#' Run RegCompass from species-aware defaults
#'
#' Prepares a pinned Human-GEM or Mouse-GEM model when `gem` is omitted and
#' delegates to the canonical GRN-first [rc_run_regcompass()] workflow.
#'
#' The workflow first normalizes single-cell RNA globally and computes ATAC
#' TF-IDF within each cell type across conditions. It then fits one Pando model
#' per condition and cell type with `peak_cor = 0.01` by default. SuperCell2
#' metacells are constructed afterwards with `gamma = 75` by default.
#'
#' @param object A Seurat RNA+ATAC object.
#' @param outdir Persistent output directory.
#' @param pfm Motif position-frequency matrices.
#' @param genome Genome object matching the selected species and ATAC coordinates.
#' @param fragment_files Must be `FALSE` for the canonical peak-count path.
#' @param species `"human"` or `"mouse"`.
#' @param gem Optional prebuilt species GEM.
#' @param gem_version Pinned model release.
#' @param medium_scenario Medium preset identifier.
#' @param medium_scenarios Optional prebuilt medium table.
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
    ...) {
  species <- match.arg(species)
  if (is.null(gem_version)) gem_version <- if (identical(species, "human")) "2.0.0" else "1.8.0"
  if (is.null(gem)) {
    gem <- rc_prepare_gem(species = species, version = gem_version)
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
    ...
  )
}
