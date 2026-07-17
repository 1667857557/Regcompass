# Workflow stage 4c: align named medium default scales with preset semantics.
#
# Human glucose/lactate presets already encode their concentration contrast in
# `uptake_fraction`. Their default scenario multiplier must therefore be one;
# applying the legacy `low_glucose = 0.1` multiplier to the complete plasma
# preset would also suppress lactate and every other available nutrient.

.rc_medium_presets_previous_default_scale <- rc_make_medium_scenarios

rc_make_medium_scenarios <- function(
    gem,
    scenario = "compass_model_bounds",
    custom_medium = NULL,
    custom_metabolites = NULL,
    uptake_scale = c(
      permissive_all_exchange = 1,
      normal_human_plasma = 1,
      rpmi1640 = 1,
      minimal = 0.1,
      low_glucose = 1,
      low_glutamine = 0.1,
      high_lactate = 1
    ),
    condition_col = NULL,
    exchange_roles = c("exchange"),
    condition = condition_col,
    exchange_limit = 1,
    strict_preset_matching = TRUE) {
  .rc_medium_presets_previous_default_scale(
    gem = gem,
    scenario = scenario,
    custom_medium = custom_medium,
    custom_metabolites = custom_metabolites,
    uptake_scale = uptake_scale,
    condition_col = condition_col,
    exchange_roles = exchange_roles,
    condition = condition,
    exchange_limit = exchange_limit,
    strict_preset_matching = strict_preset_matching
  )
}
