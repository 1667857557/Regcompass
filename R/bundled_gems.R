.rc_bundled_gem_filename <- function(spec) {
  paste0(spec$cache_prefix, "_", spec$version, "_regcompass.rds")
}

.rc_bundled_gem_path <- function(spec) {
  system.file(
    "extdata", "gem", .rc_bundled_gem_filename(spec),
    package = "RegCompassR"
  )
}

#' List GEM models bundled with RegCompassR
#'
#' @return A data frame with species, upstream release, file size, checksum,
#'   source attribution, and citation DOI.
#' @export
rc_bundled_gem_manifest <- function() {
  path <- system.file(
    "extdata", "gem", "manifest.tsv", package = "RegCompassR"
  )
  if (!nzchar(path) || !file.exists(path)) {
    return(data.frame(
      species = character(), version = character(), file = character(),
      md5 = character(), size_bytes = numeric(), source = character(),
      citation_doi = character(), stringsAsFactors = FALSE
    ))
  }
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

.rc_load_bundled_species_gem <- function(spec) {
  path <- .rc_bundled_gem_path(spec)
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  gem <- rc_read_gem(path)
  rc_validate_species_gem(gem, spec$species)
  recorded_source <- as.character(gem$model_info$source %||% "")
  recorded_version <- as.character(gem$model_info$version %||% "")
  if (!identical(recorded_source, spec$source) ||
      !identical(recorded_version, spec$version)) {
    stop(
      "Bundled GEM provenance does not match the requested source/version.",
      call. = FALSE
    )
  }
  gem$model_info$distribution <- "bundled_with_RegCompassR"
  gem$model_info$bundled_file <- basename(path)
  gem
}
