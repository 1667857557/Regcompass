\name{rc_make_medium_scenarios}
\alias{rc_make_medium_scenarios}
\title{Construct exchange-bound scenarios for a RegCompass GEM}
\usage{
rc_make_medium_scenarios(
  gem,
  scenario = "compass_model_bounds",
  custom_medium = NULL,
  uptake_scale = c(permissive_all_exchange = 1, blood_like = 1,
    culture_like = 1, minimal = 0.1, tumor_low_glucose = 0.5,
    low_glucose = 0.1, low_glutamine = 0.1, lactate_available = 1),
  condition_col = NULL,
  exchange_roles = c("exchange"),
  condition = condition_col,
  exchange_limit = 1
)
}
\arguments{
  \item{gem}{A validated RegCompass GEM, normally returned by
  \code{rc_prepare_human2_gem()}.}
  \item{scenario}{One or more scenario identifiers. The canonical default is
  \code{"compass_model_bounds"}; \code{"custom"} uses
  \code{custom_medium}.}
  \item{custom_medium}{Optional data frame with at least
  \code{medium_scenario_id}, \code{exchange_reaction_id}, \code{lb},
  \code{ub}, and \code{available}.}
  \item{uptake_scale}{A non-negative global scale or named scales for built-in
  sensitivity scenarios. It is not used to infer measured uptake rates.}
  \item{condition_col}{Deprecated compatibility alias used as the default for
  \code{condition}.}
  \item{exchange_roles}{Reaction-role labels treated as exchanges.}
  \item{condition}{Condition label for scenario rows. The canonical
  model-bound scenario must be shared and therefore uses \code{"all"}.}
  \item{exchange_limit}{Positive uniform absolute cap applied to exchange
  bounds in \code{"compass_model_bounds"}.}
}
\value{
A data frame of scenario IDs, exchange reaction IDs, lower and upper bounds,
availability flags, conditions, and bound provenance.
}
\description{
Builds validated exchange-bound tables that connect a GEM to the shared medium
used by directional RegCompass scoring.
}
\details{
The default scenario preserves the GEM's exchange directionality while capping
exchange bounds. It describes a shared model-defined environment, not measured
metabolite concentrations. Reaction IDs must match columns of \code{gem$S}.
}
\examples{
\dontrun{
gem <- rc_prepare_human2_gem("2.0.0")
medium <- rc_make_medium_scenarios(
  gem,
  scenario = "compass_model_bounds",
  exchange_limit = 1
)
head(medium)
}
}
