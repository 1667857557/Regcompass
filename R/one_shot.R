#' Run RegCompass from packaged setup defaults
#'
#' Convenience entry point for tutorials and first-pass analyses. When `gem` is
#' not supplied, the requested Human-GEM release is downloaded and converted with
#' `rc_prepare_human2_gem()`. When `medium_scenarios` is not supplied, one shared
#' medium table is built with `rc_make_medium_scenarios()` before delegating to
#' `rc_run_regcompass()`.
#' @export
rc_run_regcompass_one_shot <- function(object, outdir, pfm, genome,
                                       fragment_files = NULL,
                                       gem = NULL,
                                       humangem_version = "2.0.0",
                                       medium_scenario = "blood_like",
                                       medium_scenarios = NULL,
                                       ...) {
  if (is.null(gem)) {
    gem <- rc_prepare_human2_gem(version = humangem_version)
  }
  if (is.null(medium_scenarios)) {
    medium_scenarios <- rc_make_medium_scenarios(
      gem = gem,
      scenario = medium_scenario
    )
  }
  rc_run_regcompass(
    object = object,
    gem = gem,
    outdir = outdir,
    pfm = pfm,
    genome = genome,
    fragment_files = fragment_files,
    medium_scenarios = medium_scenarios,
    ...
  )
}
