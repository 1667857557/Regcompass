\name{rc_make_medium_scenarios}
\alias{rc_make_medium_scenarios}
\title{Construct published human and custom exchange-bound scenarios}
\usage{
rc_make_medium_scenarios(
  gem,
  scenario = "compass_model_bounds",
  custom_medium = NULL,
  custom_metabolites = NULL,
  uptake_scale = 1,
  condition_col = NULL,
  exchange_roles = c("exchange"),
  condition = condition_col,
  exchange_limit = 1,
  strict_preset_matching = TRUE
)
}
\arguments{
  \item{gem}{A validated RegCompass GEM, normally returned by
  \code{rc_prepare_human2_gem()}. Preset matching uses exchange-reaction
  annotations in \code{gem$reaction_meta}.}

  \item{scenario}{One or more scenario identifiers. Supported human presets are
  \code{"normal_human_plasma"}, \code{"high_glucose"},
  \code{"low_glucose"}, \code{"high_lactate"},
  \code{"low_lactate"}, and \code{"rpmi1640"}. The technical default is
  \code{"compass_model_bounds"}. Use \code{"custom"} with either
  \code{custom_medium} or \code{custom_metabolites}.}

  \item{custom_medium}{Optional complete reaction-level medium data frame with
  at least \code{medium_scenario_id}, \code{exchange_reaction_id},
  \code{lb}, \code{ub}, and \code{available}. Bounds are used directly and
  therefore must already be expressed in the GEM flux convention.}

  \item{custom_metabolites}{Optional metabolite-level data frame used with
  \code{scenario = "custom"}. Required columns are
  \code{metabolite_name}, \code{metabolite_pattern}, and \code{available}.
  Optional columns are \code{concentration_mM}, \code{uptake_fraction},
  \code{target_exchange_flag}, \code{required_match},
  \code{reference_label}, \code{reference_doi}, and
  \code{reference_pmid}. Patterns are matched against exchange-reaction
  annotations.}

  \item{uptake_scale}{One non-negative multiplier, or a named vector indexed by
  scenario ID, applied to preset uptake fractions. It is a sensitivity
  multiplier and is not a measured uptake rate.}

  \item{condition_col}{Deprecated compatibility alias used as the default for
  \code{condition}.}

  \item{exchange_roles}{Reaction-role labels treated as exchanges.}

  \item{condition}{Condition label written to scenario rows. Use one shared
  condition, normally \code{"all"}, when comparing biological groups unless
  measured condition-specific media are available.}

  \item{exchange_limit}{Positive absolute exchange-flux cap. For a listed
  metabolite, the uptake lower bound is limited by
  \code{-exchange_limit * uptake_fraction}; secretion is capped separately by
  the existing model and medium application logic.}

  \item{strict_preset_matching}{Logical. If \code{TRUE}, stop when required
  marker metabolites such as glucose, lactate, or glutamine cannot be matched
  to GEM exchange reactions. If \code{FALSE}, return the matched subset and
  retain matching diagnostics.}
}
\value{
A data frame of exchange-reaction bounds and provenance. Human preset rows also
contain \code{preset_metabolite}, \code{concentration_mM},
\code{uptake_fraction}, \code{reference_label}, \code{reference_doi},
\code{reference_pmid}, \code{species}, and \code{evidence_scope}. Attributes
\code{preset_diagnostics} and \code{preset_references} are attached to an
individual preset before multiple scenario tables are combined.
}
\description{
Builds exchange-bound tables that primarily determine which extracellular
metabolites may be taken up. Published presets use human plasma, human tumor or
human-cell experiments, or serum-free RPMI-1640. Users may instead provide a
metabolite availability table or exact exchange-reaction bounds.
}
\details{
When the returned table is applied with \code{rc_apply_medium_constraints()},
all exchange reactions are first closed for uptake. Only exchange reactions
listed as available in the medium table are reopened for uptake. Exchange
reactions that are absent from a preset therefore retain \code{lb = 0}.
Secretion remains allowed according to the GEM direction and
\code{rc_apply_medium_constraints()} settings.

Published concentrations are not fluxes. Values in mM are retained as
provenance. Only the designated glucose or lactate contrast is converted to a
dimensionless relative uptake cap. With the default \code{exchange_limit = 1},
the preset fractions are:

\tabular{lll}{
Preset \tab Published concentration \tab Target uptake fraction\cr
\code{normal_human_plasma} \tab glucose 5 mM; lactate 1.5 mM \tab glucose 0.20; lactate 0.075\cr
\code{high_glucose} \tab glucose 25 mM \tab glucose 1.00\cr
\code{low_glucose} \tab glucose 1 mM \tab glucose 0.04\cr
\code{high_lactate} \tab lactate 20 mM \tab lactate 1.00\cr
\code{low_lactate} \tab lactate 0.5 mM \tab lactate 0.025\cr
\code{rpmi1640} \tab glucose 11.1 mM; glutamine 2.055 mM \tab glucose 0.444; other listed nutrients available\cr
}

All other listed preset nutrients use binary availability with an uptake
fraction of 1. These values create comparable sensitivity scenarios; they do
not identify uptake in mmol per gram dry weight per hour.

The human-only literature sources encoded in each output row are:
\itemize{
  \item Normal human plasma: Cantor et al. (2017), Cell,
  DOI 10.1016/j.cell.2017.03.023, PMID 28388410; and Psychogios et al.
  (2011), PLoS One, DOI 10.1371/journal.pone.0016957, PMID 21359215.
  \item High and low glucose: Han et al. (2015), Gynecologic Oncology,
  DOI 10.1016/j.ygyno.2015.06.036, PMID 26135947; human endometrial cancer
  cells studied at 25, 5, and 1 mM glucose.
  \item High lactate: Schwickert et al. (1996), Experientia,
  DOI 10.1007/BF01919316, PMID 8641383; and Kennedy et al. (2013),
  PLoS One, DOI 10.1371/journal.pone.0075154, PMID 24069390.
  \item Low lactate: Kennedy et al. (2013), which describes a human
  physiologic blood range of approximately 0.5--2 mM.
  \item RPMI-1640: Moore et al. (1967), JAMA,
  DOI 10.1001/jama.1967.03120080053007, describing the medium for normal
  human leukocyte culture; Cantor et al. (2017) provides the human-plasma
  comparison. The preset represents serum-free basal RPMI-1640 and does not
  model fetal bovine serum.
}

Legacy names are mapped with warnings:
\code{"blood_like"} to \code{"normal_human_plasma"},
\code{"culture_like"} to \code{"rpmi1640"},
\code{"tumor_low_glucose"} to \code{"low_glucose"}, and
\code{"lactate_available"} to \code{"high_lactate"}.
}
\examples{
\dontrun{
gem <- rc_prepare_human2_gem("2.0.0")

plasma <- rc_make_medium_scenarios(
  gem,
  scenario = "normal_human_plasma",
  exchange_limit = 1
)

stress <- rc_make_medium_scenarios(
  gem,
  scenario = c("low_glucose", "high_lactate"),
  exchange_limit = 1
)

custom_metabolites <- data.frame(
  metabolite_name = c("glucose", "lactate"),
  metabolite_pattern = c("glucose|glc", "lactate|lactic acid"),
  available = TRUE,
  concentration_mM = c(3, 8),
  uptake_fraction = c(0.12, 0.40),
  target_exchange_flag = TRUE,
  required_match = TRUE,
  reference_doi = "project-specific reference"
)
custom <- rc_make_medium_scenarios(
  gem,
  scenario = "custom",
  custom_metabolites = custom_metabolites,
  exchange_limit = 1
)

# Exact measured or otherwise justified reaction-level bounds remain supported.
custom_flux <- rc_make_medium_scenarios(
  gem,
  scenario = "custom",
  custom_medium = data.frame(
    medium_scenario_id = "measured_medium",
    exchange_reaction_id = c("EXAMPLE_EXCHANGE"),
    lb = -0.2,
    ub = 1,
    available = TRUE
  )
)
}
}
\seealso{
\code{\link{rc_run_regcompass}},
\code{\link{rc_run_regcompass_one_shot}},
\code{\link{rc_prepare_human2_gem}}
}
