# Published extracellular medium scenarios

`rc_make_medium_scenarios()` accepts only medium definitions with explicit
published-paper provenance. A scenario identifier must refer to a specific
published formulation; broad biological labels, incomplete reconstructions,
manufacturer-only formulations, synthetic nutrient challenges, and technical
GEM boundary modes are not exposed as built-in scenarios.

## Current built-in scenario

| `scenario` | Species | Publication | Encoded scope |
|---|---|---|---|
| `"cantor2017_hplm"` | Human-GEM only | Cantor et al., *Cell* 2017; doi:10.1016/j.cell.2017.03.023; PMID 28388410 | HPLM components with exact published concentrations and direct one-to-one GEM exchange mapping. |

```r
medium_scenarios <- rc_make_medium_scenarios(
  gem = human_gem,
  scenario = "cantor2017_hplm",
  species = "human"
)
```

The implementation keeps only components whose concentration is directly
attributed to the Cantor 2017 HPLM formulation. It does not supplement the
medium with serum-ion values from another paper and does not infer free-ion
concentrations from ambiguous salt mixtures. Such rows are omitted rather than
approximated.

Published concentrations are provenance for extracellular availability. They
are not measured transporter fluxes and are not converted to uptake rates.
`uptake_scale` is therefore fixed at `1`. Requested exchange bounds remain the
intersection of the original GEM directionality and the shared modelling
`exchange_limit`.

## DOI-cited custom environments

Other published formulations are supplied as custom tables with
`scenario = NULL`. Every row must include a non-empty `reference_label` and a
valid published-paper `reference_doi`.

### Reaction-level input

```r
published_medium <- data.frame(
  medium_scenario_id = "published_environment_2024",
  exchange_reaction_id = c("EX_glc_D_e", "EX_gln_L_e"),
  lb = c(-0.20, -0.10),
  ub = c(1, 1),
  available = TRUE,
  reference_label = "Author et al., Journal 2024",
  reference_doi = "10.xxxx/published.article",
  stringsAsFactors = FALSE
)

medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = NULL,
  species = "human",
  custom_medium = published_medium
)
```

### Metabolite-level input

```r
published_metabolites <- data.frame(
  metabolite_name = c("glucose", "glutamine"),
  available = TRUE,
  concentration_mM = c(5, 0.55),
  uptake_fraction = 1,
  target_exchange_flag = FALSE,
  required_match = TRUE,
  reference_label = "Author et al., Journal 2024",
  reference_doi = "10.xxxx/published.article",
  stringsAsFactors = FALSE
)

medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = NULL,
  species = "human",
  custom_metabolites = published_metabolites
)
```

## Removed identifiers

The following former built-in identifiers are intentionally rejected:

```text
physiologic
normal_human_plasma
mouse_plasma
rpmi1640
dmem_high_glucose
high_glucose
low_glucose
high_lactate
low_lactate
low_glutamine
minimal
compass_model_bounds
permissive_all_exchange
custom
```

Reasons:

- `physiologic` was a broad alias rather than a paper-defined formulation;
- `normal_human_plasma` combined HPLM components with ion values from another
  source;
- `mouse_plasma` implemented only three quantitative nutrients and was not a
  complete published MPM formulation;
- RPMI/DMEM values depended partly on manufacturer formulations rather than one
  fully encoded paper definition;
- single-nutrient challenge presets grafted one published concentration onto a
  different background medium;
- `minimal`, `compass_model_bounds`, and `permissive_all_exchange` were technical
  model settings, not extracellular media;
- custom data are now passed through `scenario = NULL`, so `"custom"` no longer
  masquerades as a literature preset.

MPM or TMEM from Gardner and Stuart 2024 may be added later only after the full
published supplementary formulation is encoded and validated component by
component. Partial implementations are not retained.

## Validation and diagnostics

```r
unique(medium_scenarios[, c(
  "medium_scenario_id",
  "reference_label",
  "reference_doi",
  "evidence_scope"
)])

medium_scenarios[, c(
  "exchange_reaction_id",
  "preset_metabolite",
  "concentration_mM",
  "concentration_basis",
  "component_reference_doi",
  "lb",
  "ub",
  "rate_bound_source"
)]

attr(medium_scenarios, "preset_diagnostics")
attr(medium_scenarios, "medium_policy")
```

The expected policy is:

```text
published_paper_bound_presets_only
```

## Reference

Cantor JR, Abu-Remaileh M, Kanarek N, et al. Physiologic Medium Rewires Cellular
Metabolism and Reveals Uric Acid as an Endogenous Inhibitor of UMP Synthase.
*Cell*. 2017;169:258-272.e17. doi:10.1016/j.cell.2017.03.023.