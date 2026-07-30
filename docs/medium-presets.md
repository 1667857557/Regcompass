# Medium scenarios and published evidence

`rc_make_medium_scenarios()` exposes biological environments backed by
published plasma measurements, published culture formulations, and explicit
nutrient-challenge studies. A scenario may integrate several papers when no
single paper supplies a sufficiently complete extracellular environment. The
output records background and challenge provenance separately.

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
`permissive_all_exchange` are not biological medium scenarios and are not
accepted by this public interface. The broad alias `physiologic` is also not
accepted because it hides the species-specific evidence source.

## Plasma scenarios

| `scenario` | Species | Main sources | Quantitative policy |
|---|---|---|---|
| `normal_human_plasma` | Human-GEM | Cantor et al. 2017 HPLM; Psychogios et al. 2011 adult human serum/plasma | Exact HPLM concentrations where one-to-one; representative plasma/serum values for selected ions; unsupported rows remain availability-only. |
| `mouse_plasma` | Mouse-GEM | Gardner and Stuart 2024 mouse plasma medium; Sullivan et al. 2019 murine plasma and tumour-interstitial-fluid metabolomics | Mouse glucose, lactate and glutamine have quantitative values; the wider murine catalog is availability-only when a defensible concentration is unavailable. |

```r
human_medium <- rc_make_medium_scenarios(
  gem = human_gem,
  scenario = "normal_human_plasma",
  species = "human"
)

mouse_medium <- rc_make_medium_scenarios(
  gem = mouse_gem,
  scenario = "mouse_plasma",
  species = "mouse"
)
```

Human values are never copied into `mouse_plasma`. Murine measurements vary by
model, anatomical site, diet, and sampling method; therefore unknown mouse
concentrations remain `NA` rather than being filled from human HPLM.

## Cell-culture challenge scenarios

The five challenge presets include the common nutrients of a published basal
culture environment. The named nutrient is then replaced by the concentration
reported in the challenge paper.

| Scenario | Basal nutrient background | Target override | Challenge source |
|---|---|---:|---|
| `high_glucose` | RPMI-1640/DMEM component-availability union | glucose 25 mM | Han et al. 2015 |
| `low_glucose` | RPMI-1640/DMEM component-availability union | glucose 1 mM | Han et al. 2015 |
| `high_lactate` | DMEM nutrient background | lactate 20 mM | San-Millan et al. 2020 |
| `low_lactate` | plasma-like HPLM/Plasmax nutrient background | lactate 0.5 mM | Cho et al. 2025 |
| `low_glutamine` | DMEM nutrient background | glutamine 0.5 mM | Visagie et al. 2015 Methods |

### Why the background is retained

The target nutrient is not the only extracellular substrate available to cells.
The background therefore retains mapped amino acids, vitamins, inorganic ions,
glucose or other carbon sources, oxygen, water, and other represented nutrients
from the relevant culture formulation. This prevents a glucose, lactate, or
glutamine challenge from being interpreted as a medium containing only one
metabolite.

For `high_glucose` and `low_glucose`, Han et al. used RPMI-1640 for ECC-1 cells
and DMEM for Ishikawa cells. RegCompass stores their union as component
availability; conflicting non-target concentrations are not averaged. For
`high_lactate`, the published experiment used high-glucose DMEM with added
lactate. For `low_lactate`, the source experiment used plasma-like Plasmax.
`low_glutamine` uses the Methods-defined 0.5 mM glutamine condition; the earlier
0.05 mM package value is not retained because it is inconsistent with the
explicit Methods condition.

### Composite-scenario interpretation

A challenge scenario is a literature-backed modelling construction, not a claim
that one article supplied every row. The following columns make the construction
auditable:

```r
unique(medium_scenarios[, intersect(c(
  "medium_scenario_id",
  "medium_background_id",
  "background_reference_label",
  "background_reference_doi",
  "challenge_reference_label",
  "challenge_reference_doi",
  "scenario_construction"
), colnames(medium_scenarios))])
```

The expected construction string for these presets is:

```text
published_background_plus_named_nutrient_override
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

Each scenario receives its own medium-specific structural model. Within a
scenario, all conditions and metacells use identical exchange bounds.

## User-defined medium composition

User-defined media are supported through either `scenario = "custom"` or
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

A built-in scenario vector and one custom table may be supplied together. Custom
reaction bounds are validated for finite ordered bounds and known exchange IDs.
Custom metabolite rows are mapped through the same GEM exchange-mapping logic as
built-in media.

## Interpretation rules

1. Concentration is not uptake flux. Target concentration ratios define explicit
   relative sensitivity caps, not measured transporter rates.
2. Every requested bound is intersected with the original GEM bounds; a medium
   cannot open a direction blocked by the GEM.
3. Uptake for unlisted exchanges is closed during medium application; originally
   permitted secretion may remain open.
4. Human and mouse plasma scenarios are species-restricted.
5. Challenge scenarios are Human-GEM culture models. Mouse challenge conditions
   should be supplied as user-defined compositions with appropriate evidence.
6. Comparing `high_lactate` with `low_lactate` also compares their published
   basal backgrounds. For an isolated lactate contrast, construct two custom
   scenarios with the same background and different lactate rows.

## Diagnostics

```r
attr(medium_scenarios, "preset_diagnostics")
attr(medium_scenarios, "medium_policy")
```

The expected policy is:

```text
published_plasma_or_culture_background_with_explicit_overrides
```

## References

- Cantor JR, Abu-Remaileh M, Kanarek N, et al. *Cell*. 2017. doi:10.1016/j.cell.2017.03.023.
- Psychogios N, Hau DD, Peng J, et al. *PLoS ONE*. 2011. doi:10.1371/journal.pone.0016957.
- Gardner GL, Stuart JA. *Am J Physiol Cell Physiol*. 2024. doi:10.1152/ajpcell.00452.2024.
- Sullivan MR, Danai LV, Lewis CA, et al. *eLife*. 2019. doi:10.7554/eLife.44235.
- Moore GE, Gerner RE, Franklin HA. *JAMA*. 1967. doi:10.1001/jama.1967.03120080053007.
- Dulbecco R, Freeman G. *Virology*. 1959. doi:10.1016/0042-6822(59)90063-3.
- Vande Voorde J, Ackermann T, Pfetzer N, et al. *Science Advances*. 2019. doi:10.1126/sciadv.aau7314.
- Han J, Zhang L, Guo H, et al. *Gynecologic Oncology*. 2015. doi:10.1016/j.ygyno.2015.06.036.
- San-Millan I, Julian CG, Matarazzo C, et al. *Frontiers in Oncology*. 2020. doi:10.3389/fonc.2019.01536.
- Cho E, Spielmann G, Irving BA. *Physiological Reports*. 2025. doi:10.14814/phy2.70450.
- Visagie MH, Mqoco TV, Liebenberg L, et al. *Cell & Bioscience*. 2015. doi:10.1186/s13578-015-0030-1.