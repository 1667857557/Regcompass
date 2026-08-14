# Medium scenarios and published evidence

`rc_make_medium_scenarios()` separates biological composition metadata from metabolic flux constraints.

1. **Basal nutrient composition** comes from a high-authority reproducible formulation or quantitative extracellular metabolomics study.
2. **Challenge concentration** records the extracellular concentration used in the named glucose, lactate, or glutamine experiment.
3. **Flux bounds** are changed only by an explicit flux/bound assumption. A concentration in mM is not automatically converted into an uptake flux.

## Supported identifiers

```text
normal_human_plasma
mouse_plasma
high_glucose
low_glucose
high_lactate
low_lactate
low_glutamine
custom
```

Technical or ambiguous aliases are not biological medium presets and are not part of the public scenario API.

## Human physiological background

`normal_human_plasma` uses Cantor et al. *Cell* 2017 HPLM as the primary formulation, Rossiter et al. *Cell Metabolism* 2021 for the updated HPLM components, and Vande Voorde et al. *Science Advances* 2019 Plasmax as independent physiological-medium validation. Numerical concentrations are not averaged across publications.

```r
human_medium <- rc_make_medium_scenarios(
  gem = human_gem,
  scenario = "normal_human_plasma",
  species = "human"
)
```

The output includes the basal-medium provenance and

```text
medium_background_id = authoritative_HPLM_2017_2021
scenario_construction = authoritative_HPLM_composition_without_cross_study_averaging
```

## Mouse plasma

`mouse_plasma` is anchored to Abbott et al. *Nature* 2026. The built-in catalog is conservative: unsupported components are omitted rather than inherited from human HPLM. Glucose, lactate, and glutamine retain the secondary quantitative mouse-plasma values from Gardner and Stuart 2024; other retained components can be availability-only.

```r
mouse_medium <- rc_make_medium_scenarios(
  gem = mouse_gem,
  scenario = "mouse_plasma",
  species = "mouse"
)
```

Human concentrations are never copied into `mouse_plasma`.

## Cell-culture challenge scenarios

The five challenge presets use the same HPLM basal composition and replace only the **concentration metadata** of the named nutrient:

| Scenario | Basal background | Target concentration | Challenge source |
|---|---|---:|---|
| `high_glucose` | HPLM 2017/2021 | glucose 25 mM | Han et al. 2015 |
| `low_glucose` | HPLM 2017/2021 | glucose 1 mM | Han et al. 2015 |
| `high_lactate` | HPLM 2017/2021 | lactate 20 mM | San-Millan et al. 2020 |
| `low_lactate` | HPLM 2017/2021 | lactate 0.5 mM | Cho et al. 2025 |
| `low_glutamine` | HPLM 2017/2021 | glutamine 0.5 mM | Visagie et al. 2015 |

The canonical construction is recorded as

```text
authoritative_HPLM_background_plus_named_nutrient_concentration_metadata;no_automatic_concentration_to_flux_mapping
```

For a challenge target, the output explicitly records that concentration was not used to create the rate bound:

```text
concentration_used_for_rate_bound = FALSE
rate_bound_source = background_model_cap_concentration_metadata_only
assumption_level = literature_concentration_metadata_without_flux_conversion
```

This distinction is required dimensionally. Extracellular concentration has units of amount/volume, whereas a GEM exchange flux has units of amount/(biomass·time) or the model's equivalent normalization. Without an independently specified transporter/uptake model, there is no unique map

\[
C_{extracellular}\longrightarrow v_{uptake}.
\]

Therefore RegCompass does **not** assume `v_uptake ∝ concentration` for built-in challenge presets.

## Consequence for high/low challenge comparisons

A high- and low-concentration preset can resolve to identical final LP bounds. If so, the metabolic optimization problem is identical and the concentration metadata alone cannot produce different COMPASS/CORDA2 scores. This is intentional: the software does not manufacture a flux effect unsupported by a kinetic or experimental flux assumption.

The full-GEM cache records the actual resolved feasible-region identity using the base GEM plus canonical final bounds

\[
\{(reaction\_id,lb_{final},ub_{final})\}.
\]

Inspect:

```r
attr(full_gem_cache, "summary")[, c(
  "medium_scenario",
  "resolved_medium_fingerprint",
  "n_changed_bounds_vs_reference",
  "resolved_bounds_identical_to_reference"
)]
```

Different scenario descriptions that resolve to the same final bounds have the same `resolved_medium_fingerprint`. The input-medium fingerprint is retained only as provenance.

## Multiple scenarios

```r
medium_scenarios <- rc_make_medium_scenarios(
  gem = human_gem,
  scenario = c(
    "normal_human_plasma",
    "high_glucose",
    "low_glucose",
    "high_lactate",
    "low_lactate",
    "low_glutamine"
  ),
  species = "human"
)
```

Scenario labels and concentration metadata are retained separately even when two scenarios resolve to the same final flux bounds.

## User-defined flux assumptions

If a study supplies measured uptake limits or the analysis deliberately adopts a flux sensitivity assumption, provide that explicitly rather than deriving it implicitly from concentration.

### Exact reaction-level bounds

```r
custom_medium <- data.frame(
  medium_scenario_id = "my_measured_medium",
  exchange_reaction_id = c("EX_glc_D_e", "EX_gln_L_e"),
  lb = c(-0.20, -0.10),
  ub = c(1, 1),
  available = TRUE,
  stringsAsFactors = FALSE
)

medium_scenarios <- rc_make_medium_scenarios(
  gem = human_gem,
  scenario = "custom",
  species = "human",
  custom_medium = custom_medium
)
```

### Explicit relative uptake assumption

```r
custom_metabolites <- data.frame(
  metabolite_name = c("glucose", "glutamine", "lactate"),
  available = TRUE,
  concentration_mM = c(5, 0.55, 1.6),
  uptake_fraction = c(0.2, 0.275, 0.08),
  target_exchange_flag = TRUE,
  required_match = TRUE,
  stringsAsFactors = FALSE
)
```

Here `uptake_fraction` is an explicit modeling assumption supplied by the user. It must not be interpreted as a quantity automatically inferred from the `concentration_mM` column.

## Bound semantics

Every requested medium bound is intersected with the parent GEM bounds. A medium can restrict an exchange direction but cannot open a direction that the parent GEM blocks. Reversible exchanges are constrained through their existing lower/upper bounds; they do not need to be physically duplicated into separate reactions.

## Provenance

Inspect preset evidence with:

```r
attr(medium_scenarios, "preset_diagnostics")
attr(medium_scenarios, "medium_policy")
```

The biological-preset policy is

```text
authoritative_journal_composition_with_explicit_overrides
```

## References

- Cantor JR, Abu-Remaileh M, Kanarek N, et al. *Cell*. 2017. doi:10.1016/j.cell.2017.03.023.
- Rossiter NJ, Huggler KS, Adelmann CH, et al. *Cell Metabolism*. 2021. doi:10.1016/j.cmet.2021.02.005.
- Vande Voorde J, Ackermann T, Pfetzer N, et al. *Science Advances*. 2019. doi:10.1126/sciadv.aau7314.
- Abbott KL, Subudhi S, Ferreira R, et al. *Nature*. 2026. doi:10.1038/s41586-025-09898-9.
- Gardner GL, Stuart JA. *Am J Physiol Cell Physiol*. 2024. doi:10.1152/ajpcell.00452.2024.
- Han J, Zhang L, Guo H, et al. *Gynecologic Oncology*. 2015. doi:10.1016/j.ygyno.2015.06.036.
- San-Millan I, Julian CG, Matarazzo C, et al. *Frontiers in Oncology*. 2020. doi:10.3389/fonc.2019.01536.
- Cho E, Spielmann G, Irving BA. *Physiological Reports*. 2025. doi:10.14814/phy2.70450.
- Visagie MH, Mqoco TV, Liebenberg L, et al. *Cell & Bioscience*. 2015. doi:10.1186/s13578-015-0030-1.
