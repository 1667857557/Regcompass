# Internal compatibility aliases for retired GEM preparation helpers.
# These names remain available inside the namespace so existing internal code
# and maintenance workflows do not lose functionality. They are intentionally
# not exported and are not part of the supported user-facing API.

rc_prepare_human2_gem <- function(
    version = "2.0.0",
    cache_dir = tools::R_user_dir("RegCompassR", "cache"),
    save_rds = file.path(
      cache_dir,
      paste0("Human2_", version, "_regcompass.rds")
    ),
    force_download = FALSE,
    allow_latest = FALSE,
    source = c("auto", "bundled", "download")) {
  source <- match.arg(source)
  rc_prepare_gem(
    species = "human",
    version = version,
    cache_dir = cache_dir,
    save_rds = save_rds,
    force_download = force_download,
    allow_latest = allow_latest,
    source = source
  )
}

rc_prepare_mouse_gem <- function(
    version = "1.8.0",
    cache_dir = tools::R_user_dir("RegCompassR", "cache"),
    save_rds = file.path(
      cache_dir,
      paste0("Mouse_", version, "_regcompass.rds")
    ),
    force_download = FALSE,
    allow_latest = FALSE,
    source = c("auto", "bundled", "download")) {
  source <- match.arg(source)
  rc_prepare_gem(
    species = "mouse",
    version = version,
    cache_dir = cache_dir,
    save_rds = save_rds,
    force_download = force_download,
    allow_latest = allow_latest,
    source = source
  )
}

rc_download_species_gem <- function(...) .rc_download_species_gem(...)
rc_bundled_gem_manifest <- function() .rc_bundled_gem_manifest()
