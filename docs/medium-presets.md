# Medium scenarios and published evidence

`rc_make_medium_scenarios()` separates two evidence layers:

1. **basal nutrient composition**, which must come from a high-authority,
   reproducible formulation or quantitative extracellular metabolomics study;
2. **challenge concentration**, which may come from the experiment that defined
   the glucose, lactate, or glutamine treatment.

A challenge article is therefore not used to invent the rest of the medium. The
output records basal-composition and challenge provenance separately.

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

The technical constructions `minimal`, `compass_model_bounds`, and
`permissive_all_exchange` are not biological media. The ambiguous alias
`physiologic` is also rejected.

## Evidence hierarchy

For nutrient **composition**, the canonical priority is:

1. *Cell* or *Cell Metabolism* complete physiological formulations;
2. *Nature* quantitative plasma/interstitial-fluid metabolomics;
3. *Science Advances* independently developed physiological formulations;
4. lower-tier or older studies only for a named treatment concentration or a
   quantitatively unsupported secondary value.

No concentrations are averaged across publications. A secondary publication can
validate a formulation without contributing a numerical row.

## Human plasma and physiological culture background

`normal_human_plasma` uses:

- Cantor et al., *Cell* 2017 HPLM as the primary formulation;
- Rossiter et al., *Cell Metabolism* 2021 for the updated HPLM formulation,
  including alpha-ketoglutarate, acetylcarnitine, malate, and uridine;
- Vande Voorde et al., *Science Advances* 2019 Plasmax as independent validation
  that physiological media alter cancer-cell metabolism.

The encoded concentrations come from HPLM. Plasmax is **not numerically averaged**
with HPLM. Rounded ion values from lower-tier serum surveys are not used to fill
ambiguous salt-to-free-ion conversions.

```r
human_medium <- rc_make_medium_scenarios(
  gem = human_gem,
  scenario = "normal_human_plasma",
  species = "human"
)
```

The output includes:

```text
medium_background_id = authoritative_HPLM_2017_2021
composition_primary_reference_doi
composition_validation_reference_doi
scenario_construction = authoritative_HPLM_composition_without_cross_study_averaging
```

## Mouse plasma

`mouse_plasma` is now anchored to Abbott et al., *Nature* 2026, which quantified
absolute levels of 124 metabolites across mouse plasma, cerebrospinal fluid, and
multiple tissue interstitial fluids in NSG and C57BL/6J mice.

The built-in catalog is deliberately conservative. It retains an auditable set
of metabolites supported by that study. Components outside the supported set are
omitted rather than inherited from human HPLM. Glucose, lactate, and glutamine
retain the published mouse quantitative values from Gardner and Stuart 2024 as a
secondary source; other retained rows are availability-only unless an exact
mouse concentration is encoded.

```r
mouse_medium <- rc_make_medium_scenarios(
  gem = mouse_gem,
  scenario = "mouse_plasma",
  species = "mouse"
)
```

The output records:

```text
medium_background_id = Abbott_2026_Nature_mouse_plasma
composition_primary_reference_doi = 10.1038/s41586-025-09898-9
quantitative_secondary_reference_doi = 10.1152/ajpcell.00452.2024
```

Human concentrations are never copied into `mouse_plasma`.

## Cell-culture challenge scenarios

All five challenge presets use the **same authoritative physiological basal
composition**:

```text
Cell 2017 HPLM
+ Cell Metabolism 2021 updated HPLM components
```

Plasmax from *Science Advances* 2019 is stored as an independent validation
reference. The challenge paper overrides only the named nutrient.

| Scenario | Authoritative basal background | Target override | Challenge source |
|---|---|---:|---|
| `high_glucose` | HPLM 2017/2021 | glucose 25 mM | Han et al. 2015 |
| `low_glucose` | HPLM 2017/2021 | glucose 1 mM | Han et al. 2015 |
| `high_lactate` | HPLM 2017/2021 | lactate 20 mM | San-Millan et al. 2020 |
| `low_lactate` | HPLM 2017/2021 | lactate 0.5 mM | Cho et al. 2025 |
| `low_glutamine` | HPLM 2017/2021 | glutamine 0.5 mM | Visagie et al. 2015 Methods |

The previous RPMI/DMEM union is not used as the canonical composition source.
The original JAMA and Virology formulation papers may explain the historical
media lineage, but they are not the primary nutrient-composition evidence for
these RegCompass scenarios.

The expected construction string is:

```text
authoritative_HPLM_background_plus_named_nutrient_override
```

Inspect provenance with:

```r
unique(medium_scenarios[, intersect(c(
  "medium_scenario_id",
  "medium_background_id",
  "background_reference_label",
  "background_reference_doi",
  "background_validation_reference_label",
  "background_validation_reference_doi",
  "challenge_reference_label",
  "challenge_reference_doi",
  "scenario_construction"
), colnames(medium_scenarios))])
```

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

Each scenario receives its own medium-specific structural model. Within one
scenario, all conditions and metacells use identical exchange bounds.

## User-defined medium composition

User-defined media remain supported through either `scenario = "custom"` or
`scenario = NULL`. Publication metadata are optional and retained when supplied.

### Exact reaction-level bounds

```r
custom_medium <- data.frame(
  medium_scenario_id = "my_measured_medium",
  exchange_reaction_id = c("EX_glc_D_e", "EX_gln_L_e"),
  lb = c(-0.20, -0.10),
  ub = c(1, 1),
  available = TRUE,
  reference_label = "Optional experiment or publication label",
  reference_doi = "10.xxxx/optional.reference",
  stringsAsFactors = FALSE
)

medium_scenarios <- rc_make_medium_scenarios(
  gem = human_gem,
  scenario = "custom",
  species = "human",
  custom_medium = custom_medium
)
```

### Metabolite-level composition

```r
custom_metabolites <- data.frame(
  metabolite_name = c("glucose", "glutamine", "lactate"),
  available = TRUE,
  concentration_mM = c(5, 0.55, 1.6),
  uptake_fraction = c(0.2, 0.275, 0.08),
  target_exchange_flag = TRUE,
  required_match = TRUE,
  reference_label = "Optional experiment or publication label",
  reference_doi = "10.xxxx/optional.reference",
  stringsAsFactors = FALSE
)

medium_scenarios <- rc_make_medium_scenarios(
  gem = human_gem,
  scenario = NULL,
  species = "human",
  custom_metabolites = custom_metabolites
)
```

A built-in scenario vector and one custom table may be supplied together.

## Interpretation rules

1. **Concentration is not uptake flux.** Target concentration ratios define
   explicit relative sensitivity caps, not measured transporter rates.
2. Every requested bound is intersected with the original GEM bounds; a medium
   cannot open a direction blocked by the GEM.
3. Uptake for unlisted exchanges is closed during medium application; originally
   permitted secretion may remain open.
4. Human and mouse plasma scenarios are species-restricted.
5. The challenge scenarios use an identical HPLM background, so comparisons such
   as `high_lactate` versus `low_lactate` differ in the named target concentration
   rather than in unrelated basal nutrients.

## Diagnostics

```r
attr(medium_scenarios, "preset_diagnostics")
attr(medium_scenarios, "medium_policy")
```

The expected policy is:

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
