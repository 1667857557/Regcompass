\name{rc_run_regcompass}
\alias{rc_run_regcompass}
\title{Run the integrated RegCompass multiome workflow}
\usage{
rc_run_regcompass(
  object, gem, outdir, pfm, genome, fragment_files = NULL,
  sample_col = "sample_id", condition_col = "condition",
  celltype_col = "cell_type", rna_assay = "RNA", atac_assay = "ATAC",
  model_mode = c("meta_module_gem", "full_gem"), medium_scenarios = NULL,
  metacell_args = list(), layer1_args = list(), pando_args = list(),
  layer2_args = list(), upstream_workers = NULL, layer2_workers = NULL,
  parallel_backend = c("auto", "serial", "snow", "multicore"),
  strict_biological_defaults = TRUE,
  inference_unit = c("sample_celltype", "metacell")
)
}
\arguments{
  \item{object}{A Seurat multiome object containing RNA and ATAC assays and
  the grouping metadata columns.}
  \item{gem}{Validated GEM with a \code{gpr_table}.}
  \item{outdir}{Directory for status tables and serialized workflow outputs.}
  \item{pfm}{Position-frequency matrices passed to Pando.}
  \item{genome}{Genome object compatible with the ATAC assay and Pando.}
  \item{fragment_files}{Fragment input registered on the ATAC assay, a path or
  named list, or a data frame with \code{sample_id}, \code{assay}, and
  \code{fragment_file}.}
  \item{sample_col, condition_col, celltype_col}{Metadata column names defining
  strict condition by sample by cell-type strata.}
  \item{rna_assay, atac_assay}{Names of the RNA and chromatin assays.}
  \item{model_mode}{Use the reaction meta-module GEM or the full GEM.}
  \item{medium_scenarios}{Exchange-bound table returned by
  \code{rc_make_medium_scenarios()}. If \code{NULL}, shared COMPASS-style
  model bounds are constructed.}
  \item{metacell_args}{Named list passed to strict-stratum metacell
  construction, including settings such as \code{gamma} and minimum sizes.}
  \item{layer1_args}{Named list controlling expression/GPR calibration and
  local FASTCORE completion.}
  \item{pando_args}{Named list controlling stratum-specific Pando inference.}
  \item{layer2_args}{Named list controlling directional LP scoring. Integrated
  fields such as the GEM, unit, and medium cannot be overridden here.}
  \item{upstream_workers}{Worker count for independent strict strata.}
  \item{layer2_workers}{Worker count for shared-model scoring units.}
  \item{parallel_backend}{BiocParallel backend selection.}
  \item{strict_biological_defaults}{Logical; enforce conservative GPR and
  zero-preserving evidence defaults.}
  \item{inference_unit}{Scoring unit. \code{"sample_celltype"} is the
  confirmatory default; metacells are exploratory within-sample units.}
}
\value{
A RegCompass result list containing \code{upstream_status},
\code{upstream_barrier}, \code{layer1}, \code{grn_meta_modules},
\code{microcompass}, and recorded \code{params}. The result is also saved under
\code{outdir}.
}
\description{
Splits cells into strict condition by sample by cell-type strata, constructs
metacells and stratum-specific Pando evidence, performs local metabolic module
completion, merges sample-balanced Layer 1 evidence, builds one shared GEM, and
runs directional minimum evidence-discordance scoring.
}
\details{
The global stage runs only after retained strata finish and every biological
sample has an analyzable stratum. Inspect \code{upstream_barrier} before using
the results. The primary output is \code{microcompass$penalty}; lower values
mean stronger model support and do not represent measured flux.
}
\examples{
\dontrun{
result <- rc_run_regcompass(
  object = object, gem = gem, outdir = "RegCompass_result",
  pfm = motifs, genome = genome, fragment_files = fragment_files,
  medium_scenarios = medium,
  inference_unit = "sample_celltype"
)
stopifnot(isTRUE(result$upstream_barrier$passed))
}
}
\seealso{
\code{\link{rc_run_regcompass_one_shot}},
\code{\link{rc_prepare_human2_gem}},
\code{\link{rc_make_medium_scenarios}}
}
