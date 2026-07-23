# Predefined extracellular medium scenarios

`rc_make_medium_scenarios()` converts a named extracellular environment into reaction-level exchange bounds for the prepared Human-GEM or Mouse-GEM model. The canonical workflow uses one shared medium across conditions (`condition = "all"`) so condition contrasts are not confounded by different environmental constraints.

## Available scenarios

| `scenario` | Background and intended use | Species restrictions | Quantitative treatment |
|---|---|---|---|
| `"physiologic"` | Recommended in-vivo baseline. Resolves to `normal_human_plasma` for human and `mouse_plasma` for mouse. | Human or mouse | Literature-backed availability catalog. Glucose, lactate, and glutamine carry relative uptake caps; other listed compounds are availability constraints. |
| `"normal_human_plasma"` | Adult human plasma/HPLM-like polar nutrient environment. | Human GEM only | Human plasma/HPLM concentrations are retained as provenance; glucose, lactate, and glutamine are scaled relative to 25, 20, and 2 mM reference values. |
| `"mouse_plasma"` | Healthy-mouse plasma quantitative baseline with broader murine plasma/TIF availability support. | Mouse GEM only | Glucose 4.381 mM, lactate 3.088 mM, and glutamine 0.934 mM define relative sensitivity caps from Gardner and Stuart. Every other component is availability-only (`concentration_mM = NA`), and no human HPLM concentration or provenance is inherited. |
| `"rpmi1640"` | Serum-free RPMI-1640 basal formulation. Appropriate for an RPMI culture sensitivity analysis, not for representing serum supplementation. | Human or mouse | Formulation concentrations are stored as provenance. Availability is mapped to GEM exchanges; this is not a measured uptake-flux model. |
| `"dmem_high_glucose"` | Serum-free high-glucose DMEM basal formulation, including 25 mM glucose and 4 mM glutamine. | Human or mouse | Formulation concentrations are stored as provenance; mapped exchanges are constrained without opening directions blocked by the GEM. |
| `"high_glucose"` | Plasma background with an explicit 25 mM glucose sensitivity cap. | Uses the species-specific plasma background | Glucose uptake fraction = 25/25 = 1. Other plasma components remain unchanged. |
| `"low_glucose"` | Plasma background with an explicit 1 mM glucose sensitivity cap. | Uses the species-specific plasma background | Glucose uptake fraction = 1/25 = 0.04. |
| `"high_lactate"` | Plasma background with an explicit 20 mM lactate sensitivity cap. | Uses the species-specific plasma background | Lactate uptake fraction = 20/20 = 1. |
| `"low_lactate"` | Plasma background with an explicit 0.5 mM lactate sensitivity cap. | Uses the species-specific plasma background | Lactate uptake fraction = 0.5/20 = 0.025. |
| `"low_glutamine"` | Plasma background with an explicit 0.05 mM glutamine sensitivity cap. | Uses the species-specific plasma background | Glutamine uptake fraction = 0.05/2 = 0.025. |
| `"minimal"` | Technical minimal nutrient catalog containing glucose, glutamine, essential amino acids, oxygen, water, phosphate, bicarbonate, sodium, potassium, and chloride. | Human or mouse | Availability-only baseline with the default named `uptake_scale` of 0.1. Use as a structural sensitivity analysis, not as a physiological medium. |
| `"compass_model_bounds"` | GEM-defined environment with every annotated exchange retained in its original direction and capped uniformly. | Human or mouse | Intersects original exchange bounds with `[-exchange_limit, exchange_limit]`. |
| `"permissive_all_exchange"` | Technical all-exchange sensitivity baseline. | Human or mouse | Currently uses the same bound construction as `compass_model_bounds`, but records a different technical-assumption label. It is not a biological culture medium. |
| `"custom"` | User-defined reaction-level bounds or metabolite availability. | Must match the selected GEM | Supply exactly one of `custom_medium` or `custom_metabolites`. |

## Mouse physiological baseline

| Nutrient | Mouse concentration (mM) | High reference (mM) | `uptake_fraction` |
|---|---:|---:|---:|
| glucose | 4.381 | 25 | 0.17524 |
| lactate | 3.088 | 20 | 0.15440 |
| glutamine | 0.934 | 2 | 0.46700 |

These values come from the healthy-mouse plasma medium described by Gardner and Stuart. Sullivan et al.'s quantitative murine plasma and tumor-interstitial-fluid study supports the broader availability catalog and demonstrates that extracellular metabolite levels vary with tumor model, anatomical location, diet, and sampling. RegCompass therefore does not assign one tumor-bearing concentration to every mouse row.

## Important interpretation rules

1. **Concentration is not uptake flux.** Concentrations are retained as provenance. Only the explicitly flagged glucose, lactate, and glutamine rows convert concentration ratios into relative uptake caps; these are sensitivity assumptions, not measured transport rates.
2. **Human values do not populate the mouse preset.** Non-target mouse components are availability-only and cannot retain a human HPLM concentration, concentration basis, or component DOI.
3. **GEM directionality is never expanded.** Requested medium bounds are intersected with the original GEM bounds. A preset cannot open an uptake or secretion direction that the model originally blocked.
4. **Unlisted uptake is closed during medium application.** By default, exchanges absent from the selected catalog receive `exchange_default_lb = 0`; originally permitted secretion can remain open when `allow_secretion = TRUE`.
5. **Culture presets are basal formulations.** `rpmi1640` and `dmem_high_glucose` do not automatically represent serum, dialyzed serum, pyruvate supplementation, or laboratory-specific additives. Use `custom` for the actual experimental formulation.
6. **Compare media with all other settings fixed.** Reuse the same GEM, Layer 1 reaction expression, target reactions, target directions, solver, and convergence criteria.

## Examples

### Recommended physiological baseline

```r
medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "physiologic",
  species = "human"
)
```

### Culture-medium sensitivity analysis

```r
medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = c("rpmi1640", "dmem_high_glucose"),
  species = "human"
)
```

### Targeted nutrient sensitivity analysis

```r
medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = c("physiologic", "low_glucose", "high_lactate"),
  species = "human"
)
```

### Exact reaction-level custom medium

```r
custom_medium <- data.frame(
  medium_scenario_id = "my_medium",
  exchange_reaction_id = c("EX_glc_D_e", "EX_gln_L_e"),
  lb = c(-0.2, -0.1),
  ub = c(1, 1),
  available = TRUE
)

medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "custom",
  species = "human",
  custom_medium = custom_medium
)
```

Inspect the generated table and mapping diagnostics before scoring:

```r
unique(medium_scenarios$medium_scenario_id)
medium_scenarios[, c(
  "medium_scenario_id", "exchange_reaction_id", "preset_metabolite",
  "lb", "ub", "concentration_mM", "rate_bound_source"
)]
attr(medium_scenarios, "preset_diagnostics")
```
