\name{rc_prepare_human2_gem}
\alias{rc_prepare_human2_gem}
\title{Prepare a Human-GEM model for RegCompass}
\usage{
rc_prepare_human2_gem(
  version = "2.0.0",
  cache_dir = tools::R_user_dir("RegCompassR", "cache"),
  save_rds = file.path(cache_dir, paste0("Human2_", version, "_regcompass.rds")),
  force_download = FALSE,
  allow_latest = FALSE
)
}
\arguments{
  \item{version}{Human-GEM release version. A pinned release such as
  \code{"2.0.0"} is recommended. Use \code{"latest"} only together with
  \code{allow_latest = TRUE}.}
  \item{cache_dir}{Directory used to cache the converted model.}
  \item{save_rds}{Path of the cached RegCompass RDS model.}
  \item{force_download}{Logical; if \code{TRUE}, remove an existing cached
  model and download the requested release again.}
  \item{allow_latest}{Logical; explicitly permits the non-reproducible
  \code{"latest"} development version.}
}
\value{
A validated GEM list containing the stoichiometric matrix, reaction bounds,
reaction metadata, parsed GPR table, metabolic genes, and model provenance.
}
\description{
Downloads an official Human-GEM archive, converts its YAML model to the
RegCompass representation, attaches gene-protein-reaction and reaction-role
annotations, validates the result, and caches it as an RDS file.
}
\details{
The returned reaction identifiers are the identifiers that must be used in
medium, target-reaction, and reaction-membership inputs. Pin \code{version} and
retain \code{model_info} for reproducible analyses.
}
\examples{
\dontrun{
gem <- rc_prepare_human2_gem(version = "2.0.0")
gem$model_info[c("version", "commit", "checksum")]
head(gem$gpr_table)
}
}
