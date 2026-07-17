# Supported functions

Only four functions form the supported API. Internal helpers may change and
should not be called directly.

After installation, each supported function has a standard R help page. Open
it in RStudio with `?rc_run_regcompass`, `?rc_run_regcompass_one_shot`,
`?rc_prepare_human2_gem`, or `?rc_make_medium_scenarios`.

## `rc_run_regcompass_one_shot()`

The tutorial entry point. It prepares a Human-GEM model and shared medium when
needed, then calls `rc_run_regcompass()`.

```r
result <- rc_run_regcompass_one_shot(
  object, "RegCompass_result", motifs, genome, fragment_files,
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type"
)
```

## `rc_prepare_human2_gem()`

Downloads and converts one Human-GEM release.

```r
gem <- rc_prepare_human2_gem(version = "2.0.0")
```

## `rc_make_medium_scenarios()`

Creates exchange constraints used by every condition. The default
`compass_model_bounds` remains a shared technical baseline that preserves GEM
directionality and caps exchange reactions. It is not a measured biological
medium.

```r
technical_medium <- rc_make_medium_scenarios(
  gem,
  scenario = "compass_model_bounds",
  exchange_limit = 1
)
```

### Published human presets

The function now provides six human-only presets:

| Preset | Main definition | Human references |
|---|---|---|
| `normal_human_plasma` | Plasma/serum metabolite availability; glucose 5 mM and lactate 1.5 mM | Cantor et al., *Cell* 2017, DOI `10.1016/j.cell.2017.03.023`; Psychogios et al., *PLoS One* 2011, DOI `10.1371/journal.pone.0016957` |
| `high_glucose` | Human plasma background with glucose 25 mM | Han et al., *Gynecologic Oncology* 2015, DOI `10.1016/j.ygyno.2015.06.036` |
| `low_glucose` | Human plasma background with glucose 1 mM | Han et al., *Gynecologic Oncology* 2015, DOI `10.1016/j.ygyno.2015.06.036` |
| `high_lactate` | Human plasma background with lactate 20 mM | Schwickert et al., *Experientia* 1996, DOI `10.1007/BF01919316`; Kennedy et al., *PLoS One* 2013, DOI `10.1371/journal.pone.0075154` |
| `low_lactate` | Human plasma background with lactate 0.5 mM | Kennedy et al., *PLoS One* 2013, DOI `10.1371/journal.pone.0075154` |
| `rpmi1640` | Serum-free basal RPMI-1640 nutrient availability; glucose 11.1 mM and glutamine 2.055 mM | Moore et al., *JAMA* 1967, DOI `10.1001/jama.1967.03120080053007`; Cantor et al., *Cell* 2017 |

```r
plasma <- rc_make_medium_scenarios(
  gem,
  scenario = "normal_human_plasma"
)

nutrient_stress <- rc_make_medium_scenarios(
  gem,
  scenario = c("low_glucose", "high_lactate")
)

culture <- rc_make_medium_scenarios(
  gem,
  scenario = "rpmi1640"
)
```

The presets primarily implement **allow/deny uptake**. When a medium table is
applied, exchange uptake is closed first; only exchange reactions represented
by available preset metabolites are reopened. Unlisted exchanges therefore
remain at `lb = 0`. Secretion remains governed by the GEM and the medium
application settings.

Concentrations are not fluxes. They are retained in `concentration_mM` as
provenance. Only the designated glucose or lactate contrast is converted to a
dimensionless relative uptake fraction. With `exchange_limit = 1`, the target
fractions are 1.00 for 25 mM glucose or 20 mM lactate, 0.20 for normal plasma
glucose, 0.075 for normal plasma lactate, 0.04 for 1 mM glucose, and 0.025 for
0.5 mM lactate. These are sensitivity bounds, not measured
`mmol / gDW / h` rates.

Each human preset output includes its paper citation, DOI, PMID, species,
concentration, matching diagnostics, and bound provenance. The function stops
by default if required marker metabolites cannot be matched to Human-GEM
exchange annotations.

### User-defined environments

Users can provide either a complete exchange-reaction table or a metabolite
availability table.

```r
# Metabolite-level input: the function maps patterns to exchange reactions.
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

# Reaction-level input: use exact measured or justified model bounds directly.
custom_flux <- rc_make_medium_scenarios(
  gem,
  scenario = "custom",
  custom_medium = data.frame(
    medium_scenario_id = "measured_medium",
    exchange_reaction_id = "EXAMPLE_EXCHANGE",
    lb = -0.2,
    ub = 1,
    available = TRUE
  )
)
```

Current named backgrounds include `normal_human_plasma`, `rpmi1640`,
`low_glucose`, and `high_lactate`. Retired names are not compatibility aliases.
Use `custom_medium` only when measured or justified bounds are available.

## `rc_run_regcompass()`

Runs the canonical workflow with an explicit GEM and medium:

```r
result <- rc_run_regcompass(
  object, gem, "RegCompass_result", motifs, genome, fragment_files,
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type",
  medium_scenarios = medium
)
```

The main path is strict-stratum metacells and Pando inference, local FASTCORE,
global calibration and a shared GEM, followed by directional scoring. The
[parameter-selection guide](../README.md#choosing-analysis-parameters) explains
how to choose metacell resolution, Pando thresholds, calibration options, and
the LP solver. Project code should handle downstream reporting and statistics.
