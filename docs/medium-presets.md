# Predefined extracellular medium scenarios

`rc_make_medium_scenarios()` converts a named extracellular environment into reaction-level exchange bounds for a prepared Human-GEM or Mouse-GEM model. Preset concentrations are provenance for extracellular availability and sensitivity analysis; they are **not measured uptake fluxes**.

## Species policy

The biological origin of a concentration and the species accepted by the software are kept explicit.

| `scenario` | Accepted GEM species | Biological or formulation origin | Intended use |
|---|---|---|---|
| `"physiologic"` | Human or mouse | Resolves to `normal_human_plasma` for Human-GEM and `mouse_plasma` for Mouse-GEM | Recommended in-vivo baseline |
| `"normal_human_plasma"` | **Human only** | Adult human plasma/HPLM | Human physiological extracellular environment |
| `"mouse_plasma"` | **Mouse only** | Healthy-mouse plasma medium plus murine plasma/TIF availability evidence | Mouse physiological extracellular environment |
| `"rpmi1640"` | Human or mouse | Serum-free RPMI-1640 chemical formulation | Culture-medium sensitivity analysis; not species-specific physiology |
| `"dmem_high_glucose"` | Human or mouse | Serum-free high-glucose DMEM chemical formulation | Culture-medium sensitivity analysis; not species-specific physiology |
| `"high_glucose"` | **Human only** | Human-cell glucose challenge at 25 mM | High-glucose sensitivity on the human-plasma background |
| `"low_glucose"` | **Human only** | Human-cell glucose challenge at 1 mM | Low-glucose sensitivity on the human-plasma background |
| `"high_lactate"` | **Human only** | Human MCF7-cell lactate challenge at 20 mM | High-lactate sensitivity on the human-plasma background |
| `"low_lactate"` | **Human only** | Human T-cell lactate condition at 0.5 mM | Low-lactate sensitivity on the human-plasma background |
| `"low_glutamine"` | **Human only** | Human breast/cervical tumour-cell glutamine deprivation at 0.05 mM | Low-glutamine sensitivity on the human-plasma background |
| `"minimal"` | Human or mouse | Technical nutrient catalog | Structural sensitivity analysis, not a physiological medium |
| `"compass_model_bounds"` | Human or mouse | Original GEM exchange directions with a uniform cap | Model-defined technical baseline |
| `"permissive_all_exchange"` | Human or mouse | Current all-exchange technical construction | Technical sensitivity baseline |
| `"custom"` | Must match the selected GEM | User supplied | Measured or explicitly assumed environment |

The five nutrient-challenge presets were previously allowed to inherit `mouse_plasma`. Their defining concentrations are from human experimental systems, so that behaviour mixed species evidence. They now require Human-GEM. For a mouse nutrient challenge, use `scenario = "custom"` with a murine or experiment-specific concentration.

RPMI-1640 and DMEM are not “human values” or “mouse values”: they are chemically defined culture formulations. They remain usable with either GEM, but do not represent serum supplementation, dialysed serum, pyruvate additions, or laboratory-specific modifications.

## Quantitative provenance

### Physiological presets

`normal_human_plasma` uses two evidence layers:

1. Components mapping one-to-one to Human Plasma-Like Medium (HPLM) use the exact formulation concentration. Examples include glucose 5.0 mM, lactate 1.6 mM, glutamine 0.55000347 mM, alanine 0.43000337 mM, serine 0.15000476 mM, pyruvate 0.049999997 mM, acetate 0.039997563 mM, and citrate 0.13000052 mM.
2. Ions represented by several salts in HPLM retain representative adult human serum/plasma values from Psychogios et al. rather than inferring a free-ion concentration by summing formulation salts.

The reaction-level output records `concentration_mM`, `concentration_basis`, and `component_reference_doi` for every mapped quantitative component. Components without defensible quantitative correspondence remain availability-only (`concentration_mM = NA`).

`mouse_plasma` does not inherit any human HPLM concentration. Only the three target nutrients below receive quantitative sensitivity caps; other murine rows remain availability-only because plasma and tumour-interstitial-fluid concentrations depend on tumour model, anatomical site, diet, and sampling procedure.

| Nutrient | Human HPLM baseline (mM) | Mouse plasma-medium baseline (mM) | High reference for relative cap (mM) |
|---|---:|---:|---:|
| glucose | 5.0 | 4.381 | 25 |
| lactate | 1.6 | 3.088 | 20 |
| glutamine | 0.55000347 | 0.934 | 2 |

The relative cap is `min(1, concentration / high_reference)`. It is a modelling sensitivity assumption, not a measured transporter rate.

### Culture formulations

| Preset | Defining formulation values | Source |
|---|---|---|
| `rpmi1640` | glucose 11.111111 mM; glutamine 2.0547945 mM; other stored components follow Gibco formulation 11875 | Moore et al. and the manufacturer formulation |
| `dmem_high_glucose` | glucose 25 mM; glutamine 4 mM; other stored components follow Gibco formulation 11965 | Dulbecco and Freeman and the manufacturer formulation |

These presets describe basal formulation only. Protein, lipid, hormone, and growth-factor contributions from serum are not inferred.

### Human nutrient challenges

| Preset | Overridden nutrient | Concentration (mM) | Experimental provenance | DOI |
|---|---|---:|---|---|
| `high_glucose` | glucose | 25 | ECC-1 and Ishikawa human endometrial cancer cells | [10.1016/j.ygyno.2015.06.036](https://doi.org/10.1016/j.ygyno.2015.06.036) |
| `low_glucose` | glucose | 1 | ECC-1 and Ishikawa human endometrial cancer cells | [10.1016/j.ygyno.2015.06.036](https://doi.org/10.1016/j.ygyno.2015.06.036) |
| `high_lactate` | lactate | 20 | Human MCF7 breast cancer cells | [10.3389/fonc.2019.01536](https://doi.org/10.3389/fonc.2019.01536) |
| `low_lactate` | lactate | 0.5 | Human T cells | [10.14814/phy2.70450](https://doi.org/10.14814/phy2.70450) |
| `low_glutamine` | glutamine | 0.05 | Human breast and cervical tumourigenic cell lines | [10.1186/s13578-015-0030-1](https://doi.org/10.1186/s13578-015-0030-1) |

Only the named nutrient is replaced. The remaining components come from `normal_human_plasma`. The overridden row records the challenge study in `component_reference_doi` rather than retaining the baseline plasma DOI.

## Interpretation rules

1. **Concentration is not uptake flux.** Only explicitly flagged glucose, lactate, and glutamine rows convert a concentration ratio into a relative uptake cap.
2. **Species provenance is enforced.** Human plasma and human nutrient challenges cannot be used with Mouse-GEM; mouse plasma cannot be used with Human-GEM.
3. **Culture formulations are species-neutral, not physiological.** Their values describe bottle formulations rather than an organism.
4. **GEM directionality is never expanded.** Requested bounds are intersected with original GEM bounds.
5. **Unlisted uptake is closed during medium application.** Exchanges absent from the catalog receive `exchange_default_lb = 0`; originally permitted secretion may remain open when `allow_secretion = TRUE`.
6. **Use a custom medium for the actual experiment.** This includes serum supplementation, dialysed serum, added pyruvate, mouse-specific nutrient challenges, measured plasma/TIF samples, and laboratory-specific formulations.

## Examples

```r
human_medium <- rc_make_medium_scenarios(
  gem = human_gem,
  scenario = "physiologic",
  species = "human"
)

mouse_medium <- rc_make_medium_scenarios(
  gem = mouse_gem,
  scenario = "physiologic",
  species = "mouse"
)

human_low_glucose <- rc_make_medium_scenarios(
  gem = human_gem,
  scenario = "low_glucose",
  species = "human"
)
```

For a mouse experiment-specific challenge:

```r
mouse_glucose <- data.frame(
  medium_scenario_id = "mouse_glucose_measured",
  exchange_reaction_id = "EX_glc_D_e",
  lb = -0.2,
  ub = 1,
  available = TRUE
)

mouse_medium <- rc_make_medium_scenarios(
  gem = mouse_gem,
  scenario = "custom",
  species = "mouse",
  custom_medium = mouse_glucose
)
```

Inspect provenance and mapping diagnostics before scoring:

```r
medium_scenarios[, c(
  "medium_scenario_id", "exchange_reaction_id", "preset_metabolite",
  "lb", "ub", "concentration_mM", "concentration_basis",
  "component_reference_doi", "rate_bound_source"
)]
attr(medium_scenarios, "preset_diagnostics")
```

## References

- Cantor JR, Abu-Remaileh M, Kanarek N, et al. Physiologic Medium Rewires Cellular Metabolism and Reveals Uric Acid as an Endogenous Inhibitor of UMP Synthase. *Cell*. 2017;169:258-272.e17. [doi:10.1016/j.cell.2017.03.023](https://doi.org/10.1016/j.cell.2017.03.023).
- Psychogios N, Hau DD, Peng J, et al. The Human Serum Metabolome. *PLoS ONE*. 2011;6:e16957. [doi:10.1371/journal.pone.0016957](https://doi.org/10.1371/journal.pone.0016957).
- Gardner GL, Stuart JA. Tumor microenvironment-like conditions alter pancreatic cancer cell metabolism and behavior. *Am J Physiol Cell Physiol*. 2024. [doi:10.1152/ajpcell.00452.2024](https://doi.org/10.1152/ajpcell.00452.2024).
- Sullivan MR, Danai LV, Lewis CA, et al. Quantification of microenvironmental metabolites in murine cancers reveals determinants of tumor nutrient availability. *eLife*. 2019;8:e44235. [doi:10.7554/eLife.44235](https://doi.org/10.7554/eLife.44235).
- Moore GE, Gerner RE, Franklin HA. Culture of normal human leukocytes. *JAMA*. 1967;199:519-524. [doi:10.1001/jama.1967.03120080053007](https://doi.org/10.1001/jama.1967.03120080053007).
- Dulbecco R, Freeman G. Plaque production by the polyoma virus. *Virology*. 1959;8:396-397. [doi:10.1016/0042-6822(59)90063-3](https://doi.org/10.1016/0042-6822(59)90063-3).
- Han J, Zhang L, Guo H, et al. Glucose promotes cell proliferation, glucose uptake and invasion in endometrial cancer cells via AMPK/mTOR/S6 and MAPK signaling. *Gynecol Oncol*. 2015;138:668-675. [doi:10.1016/j.ygyno.2015.06.036](https://doi.org/10.1016/j.ygyno.2015.06.036).
- San-Millán I, Julian CG, Matarazzo C, Martinez J, Brooks GA. Is Lactate an Oncometabolite? *Front Oncol*. 2019;9:1536. [doi:10.3389/fonc.2019.01536](https://doi.org/10.3389/fonc.2019.01536).
- Cho E, Spielmann G, Irving BA. Effects of lactate concentration on T-cell phenotype and mitochondrial respiration. *Physiol Rep*. 2025;13:e70450. [doi:10.14814/phy2.70450](https://doi.org/10.14814/phy2.70450).
- Visagie MH, Mqoco TV, Liebenberg L, et al. Influence of partial and complete glutamine- and glucose-deprivation of breast- and cervical tumourigenic cell lines. *Cell Biosci*. 2015;5:37. [doi:10.1186/s13578-015-0030-1](https://doi.org/10.1186/s13578-015-0030-1).
- Gibco/Thermo Fisher Scientific. [Human Plasma-Like Medium formulation](https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.360.html), [RPMI-1640 formulation 11875](https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.114.html), and [DMEM high-glucose formulation 11965](https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.8.html).
