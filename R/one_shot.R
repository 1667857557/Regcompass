#' Run RegCompass from species-aware defaults
#'
#' Downloads a pinned Human-GEM or Mouse-GEM release when `gem` is not supplied,
#' builds a species-matched literature-backed physiological medium by default,
#' and delegates to `rc_run_regcompass()`.
#'
#' @param object A Seurat RNA+ATAC object.
#' @param outdir Persistent output directory.
#' @param pfm Motif position-frequency matrices.
#' @param genome Genome object matching the selected species and ATAC coordinates.
#' @param fragment_files Optional fragment-file mapping.
#' @param species `"human"` or `"mouse"`.
#' @param gem Optional prebuilt species GEM.
#' @param gem_version Pinned model release. Defaults to Human-GEM 2.0.0 or
#'   Mouse-GEM 1.8.0.
#' @param humangem_version Deprecated Human-GEM compatibility argument.
#' @param medium_scenario Medium preset identifier. The default `"physiologic"`
#'   resolves to human or mouse plasma.
#' @param medium_scenarios Optional prebuilt medium table.
#' @param ... Arguments passed to `rc_run_regcompass()`.
#' @return A RegCompass result list.
#' @export
rc_run_regcompass_one_shot <- function(
    object, outdir, pfm, genome,
    fragment_files = NULL,
    species = c("human", "mouse"),
    gem = NULL,
    gem_version = NULL,
    humangem_version = NULL,
    medium_scenario = "physiologic",
    medium_scenarios = NULL,
    ...) {
  species <- match.arg(species)
  if (!is.null(humangem_version)) {
    if (!identical(species, "human")) {
      stop("`humangem_version` cannot be used with `species = 'mouse'`.",
           call. = FALSE)
    }
    warning(
      "`humangem_version` is deprecated; use `gem_version`.",
      call. = FALSE
    )
    if (is.null(gem_version)) gem_version <- humangem_version
  }
  if (is.null(gem_version)) {
    gem_version <- if (identical(species, "human")) "2.0.0" else "1.8.0"
  }
  if (is.null(gem)) {
    gem <- rc_prepare_gem(
      species = species,
      version = gem_version
    )
  } else {
    .rc_infer_gem_species(gem, species)
  }
  if (is.null(medium_scenarios)) {
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
    pfm = pfm,
    genome = genome,
    fragment_files = fragment_files,
    species = species,
    medium_scenarios = medium_scenarios,
    ...
  )
}
