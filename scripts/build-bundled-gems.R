options(timeout = 18000)
dir.create("inst/extdata/gem", recursive = TRUE, showWarnings = FALSE)

collate <- read.dcf("DESCRIPTION", fields = "Collate")[[1L]]
source_order <- scan(text = collate, what = character(), quiet = TRUE)
source_order <- gsub("^'|'$", "", source_order)
invisible(lapply(file.path("R", source_order), source, local = .GlobalEnv))

build_one <- function(species, version, filename) {
  temporary <- tempfile(fileext = ".rds")
  gem <- rc_prepare_gem(
    species = species,
    version = version,
    save_rds = temporary,
    force_download = TRUE,
    source = "download"
  )
  gem$model_info$distribution <- "bundled_with_RegCompassR"
  gem$model_info$bundled_package_version <- "1.8.3"
  output <- file.path("inst/extdata/gem", filename)
  saveRDS(gem, output, compress = "xz")
  rc_validate_species_gem(readRDS(output), species)
  data.frame(
    species = species,
    version = version,
    file = filename,
    md5 = unname(tools::md5sum(output)),
    size_bytes = file.info(output)$size,
    source = as.character(gem$model_info$source),
    citation_doi = as.character(gem$model_info$citation_doi),
    stringsAsFactors = FALSE
  )
}

manifest <- rbind(
  build_one("human", "2.0.0", "Human2_2.0.0_regcompass.rds"),
  build_one("mouse", "1.8.0", "Mouse_1.8.0_regcompass.rds")
)
utils::write.table(
  manifest,
  file = "inst/extdata/gem/manifest.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)
