options(timeout = 18000)
dir.create("inst/extdata/gem", recursive = TRUE, showWarnings = FALSE)

collate <- read.dcf("DESCRIPTION", fields = "Collate")[[1L]]
source_order <- scan(text = collate, what = character(), quiet = TRUE)
source_order <- gsub("^'|'$", "", source_order)
invisible(lapply(file.path("R", source_order), source, local = .GlobalEnv))

bundle_release_date <- "2026-07-24"

download_checked <- function(url, destination) {
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  status <- utils::download.file(url, destination, mode = "wb", quiet = FALSE)
  if (!identical(as.integer(status), 0L) || !file.exists(destination) ||
      is.na(file.info(destination)$size) || file.info(destination)$size <= 0) {
    stop("Failed to download required GEM source file: ", url, call. = FALSE)
  }
  destination
}

prepare_minimal_source <- function(species, version) {
  spec <- .rc_species_gem_spec(species, version)
  ref <- paste0("v", version)
  root <- tempfile(paste0(spec$repository_name, "-minimal-"))
  model_dir <- file.path(root, "model")
  base <- paste0(
    "https://raw.githubusercontent.com/SysBioChalmers/",
    spec$repository_name, "/", ref, "/model/"
  )
  download_checked(
    paste0(base, spec$model_file),
    file.path(model_dir, spec$model_file)
  )
  download_checked(
    paste0(base, "reactions.tsv"),
    file.path(model_dir, "reactions.tsv")
  )
  if (identical(species, "human")) {
    download_checked(
      paste0(base, "genes.tsv"),
      file.path(model_dir, "genes.tsv")
    )
  }
  list(root = root, spec = spec, ref = ref)
}

build_one <- function(species, version, filename) {
  source <- prepare_minimal_source(species, version)
  prepared <- rc_prepare_species_gpr_table(
    source$root,
    species = species,
    gene_format = source$spec$gene_format
  )
  model_yml <- file.path(source$root, "model", source$spec$model_file)
  checksum <- unname(tools::md5sum(model_yml)[[1L]])
  gem <- rc_convert_yaml_to_regcompass(
    model_yml = model_yml,
    species = species,
    version = version,
    commit = source$ref,
    checksum = checksum
  )
  gem <- rc_enrich_humangem_metadata(
    gem,
    reactions_tsv = prepared$reactions,
    model_yml = model_yml
  )
  gem <- rc_annotate_reaction_roles(gem, overwrite_existing = TRUE)
  gem$gpr_table <- prepared$gpr_table
  gem$metabolic_genes <- prepared$metabolic_genes
  gem$reaction_rules <- prepared$reaction_rules
  gem$genes <- prepared$genes
  gem$reactions <- prepared$reactions
  gem$model_info$conversion_date <- bundle_release_date
  gem$model_info$gene_format <- source$spec$gene_format
  gem$model_info$archive <- NA_character_
  gem$model_info$archive_url <- NA_character_
  gem$model_info$annotation_schema <- "regcompass_species_gem_v1"
  gem$model_info$citation <- source$spec$citation
  gem$model_info$citation_doi <- source$spec$citation_doi
  gem$model_info$distribution <- "bundled_with_RegCompassR"
  gem$model_info$bundled_package_version <- "1.8.3"
  gem$model_info$upstream_ref <- source$ref
  gem$model_info$source_files <- c(
    source$spec$model_file,
    "reactions.tsv",
    if (identical(species, "human")) "genes.tsv" else character()
  )
  rc_validate_species_gem(gem, species)

  output <- file.path("inst/extdata/gem", filename)
  saveRDS(gem, output, compress = "xz", version = 3)
  reloaded <- readRDS(output)
  rc_validate_species_gem(reloaded, species)
  if (!identical(.rc_stage_gem_fingerprint(reloaded),
                 .rc_stage_gem_fingerprint(gem))) {
    stop("Bundled GEM fingerprint changed after serialization.", call. = FALSE)
  }
  data.frame(
    species = species,
    version = version,
    file = filename,
    md5 = unname(tools::md5sum(output)),
    size_bytes = as.numeric(file.info(output)$size),
    source = as.character(gem$model_info$source),
    upstream_ref = source$ref,
    upstream_model_md5 = checksum,
    conversion_date = bundle_release_date,
    citation_doi = as.character(gem$model_info$citation_doi),
    license = "CC-BY-4.0",
    stringsAsFactors = FALSE
  )
}

manifest <- rbind(
  build_one("human", "2.0.0", "Human2_2.0.0_regcompass.rds"),
  build_one("mouse", "1.8.0", "Mouse_1.8.0_regcompass.rds")
)
if (any(manifest$size_bytes >= 100 * 1024^2)) {
  stop("A bundled GEM exceeds GitHub's 100 MB single-file limit.", call. = FALSE)
}
utils::write.table(
  manifest,
  file = "inst/extdata/gem/manifest.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)
