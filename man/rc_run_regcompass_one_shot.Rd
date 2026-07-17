\name{rc_run_regcompass_one_shot}
\alias{rc_run_regcompass_one_shot}
\title{Run RegCompass with packaged setup defaults}
\usage{
rc_run_regcompass_one_shot(
  object, outdir, pfm, genome, fragment_files = NULL, gem = NULL,
  humangem_version = "2.0.0", medium_scenario = "compass_model_bounds",
  medium_scenarios = NULL, ...
)
}
\arguments{
  \item{object}{A Seurat multiome object accepted by
  \code{rc_run_regcompass()}.}
  \item{outdir}{Directory for status tables and serialized outputs.}
  \item{pfm}{Position-frequency matrices passed to Pando.}
  \item{genome}{Genome object matching the ATAC assay.}
  \item{fragment_files}{Fragment inputs accepted by
  \code{rc_run_regcompass()}.}
  \item{gem}{Optional prepared GEM. If \code{NULL}, the requested Human-GEM
  release is downloaded and converted.}
  \item{humangem_version}{Pinned Human-GEM release used when \code{gem} is
  \code{NULL}.}
  \item{medium_scenario}{Scenario passed to
  \code{rc_make_medium_scenarios()} when \code{medium_scenarios} is
  \code{NULL}. In addition to \code{"compass_model_bounds"}, published
  human presets include \code{"normal_human_plasma"},
  \code{"high_glucose"}, \code{"low_glucose"},
  \code{"high_lactate"}, \code{"low_lactate"}, and
  \code{"rpmi1640"}. See \code{?rc_make_medium_scenarios} for values,
  citations and uptake semantics.}
  \item{medium_scenarios}{Optional pre-built medium table. When supplied, it
  takes precedence over \code{medium_scenario}. Build this table explicitly
  when using \code{custom_metabolites}, exact reaction bounds, a non-default
  \code{exchange_limit}, or multiple scenarios.}
  \item{...}{Additional arguments passed to \code{rc_run_regcompass()}, such as
  metadata columns, argument bundles, workers, and \code{inference_unit}.}
}
\value{
The same integrated result list returned by \code{rc_run_regcompass()}.
}
\description{
Convenience entry point that prepares a pinned Human-GEM model and shared
medium when they are not supplied, then delegates to the canonical integrated
RegCompass workflow.
}
\details{
The one-shot helper accepts one named medium preset through
\code{medium_scenario}. For full control, first call
\code{rc_make_medium_scenarios()} and pass its output through
\code{medium_scenarios}. Published concentrations are provenance and relative
sensitivity inputs, not physical uptake rates.
}
\examples{
\dontrun{
result <- rc_run_regcompass_one_shot(
  object = object,
  outdir = "RegCompass_result",
  pfm = motifs,
  genome = genome,
  fragment_files = fragment_files,
  medium_scenario = "normal_human_plasma",
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type",
  inference_unit = "sample_celltype"
)
}
}
\seealso{
\code{\link{rc_run_regcompass}},
\code{\link{rc_prepare_human2_gem}},
\code{\link{rc_make_medium_scenarios}}
}
