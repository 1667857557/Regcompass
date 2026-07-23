# Plasma medium evidence and GEM mapping

## Scope

The physiological presets represent a curated, metabolism-relevant **polar
small-molecule availability catalog**. They are not an exhaustive reconstruction
of plasma, serum, lipoproteins, protein-bound nutrients, hormones, or growth
factors.

`normal_human_plasma` combines:

- the Human Plasma-Like Medium (HPLM) formulation and its measured plasma
  concentrations;
- the human serum metabolome reported by Psychogios et al.;
- inorganic ions and gases required by the GEM.

`mouse_plasma` uses two explicitly separated evidence layers:

1. **Quantitative healthy-mouse baseline.** Gardner and Stuart's mouse plasma
   medium provides the concentrations used for the only three nutrients that
   RegCompass converts into relative uptake caps: glucose 4.381 mM, lactate
   3.088 mM, and glutamine 0.934 mM.
2. **Broader murine availability evidence.** Sullivan et al. quantified more
   than 118 metabolites in mouse plasma and tumor interstitial fluid and showed
   that local concentrations depend on tumor type, anatomical location, diet,
   and sampling context. These measurements support the broader availability
   catalog, but are not collapsed into one universal quantitative medium.

All mouse components other than glucose, lactate, and glutamine are therefore
availability-only rows (`concentration_mM = NA`). They receive no quantitative
uptake cap and never inherit a human HPLM concentration or provenance label.

## Quantitative target nutrients

| RegCompass name | Mouse reference concentration (mM) | Existing high reference (mM) | Relative uptake fraction | Evidence |
|---|---:|---:|---:|---|
| glucose | 4.381 | 25 | 0.17524 | Gardner and Stuart 2024 |
| lactate | 3.088 | 20 | 0.15440 | Gardner and Stuart 2024 |
| glutamine | 0.934 | 2 | 0.46700 | Gardner and Stuart 2024 |

The fractions are sensitivity assumptions calculated as reference concentration
divided by the existing RegCompass high-concentration reference. They are not
measured membrane-transport rates. The resulting requested exchange bounds are
still intersected with the original GEM directionality.

## Added polar nutrients

| RegCompass name | GEM exact aliases | Human reference concentration (mM) | Evidence |
|---|---|---:|---|
| glutathione | glutathione; reduced glutathione; GSH | 0.024999188 | HPLM |
| acetone | acetone | 0.06000344 | HPLM |
| fructose | D-fructose; fructose | 0.03999778 | HPLM |
| galactose | D-galactose; galactose | 0.060002223 | HPLM |
| glycerol | glycerol | 0.11999132 | HPLM |
| hypoxanthine | hypoxanthine | 0.010000000 | HPLM |
| uridine | uridine | 0.003001638 | HPLM |
| acetylcarnitine | O-acetylcarnitine; L-acetylcarnitine | 0.005002086 | HPLM |
| betaine | betaine; trimethylglycine | 0.070004270 | HPLM |
| alpha_aminobutyrate | L-2-aminobutanoate; 2-aminobutyrate | 0.019996122 | HPLM |
| citrulline | L-citrulline; citrulline | 0.040000000 | HPLM |
| ornithine | L-ornithine; ornithine | 0.069997630 | HPLM |
| n_acetylglycine | N-acetylglycine; acetylglycine | 0.089999996 | HPLM |
| taurine | taurine | 0.090003190 | HPLM |
| alpha_ketoglutarate | 2-oxoglutarate; alpha-ketoglutarate; AKG | 0.005003422 | HPLM |
| formate | formate; formic acid | 0.049989140 | HPLM |
| malate | L-malate; malate | 0.0049966443 | HPLM |
| malonate | malonate; malonic acid | 0.010003844 | HPLM |
| succinate | succinate; succinic acid | 0.020001695 | HPLM |
| ammonium | ammonium; NH4+ | 0.12002288 | HPLM formulation |
| nitrate | nitrate; NO3- | 0.01999783 | HPLM formulation, stoichiometric nitrate content |

The human concentration column is provenance for `normal_human_plasma` only.
Except for explicitly flagged glucose/lactate/glutamine sensitivity rows,
RegCompass does not interpret a medium concentration as a measured uptake flux.

## Exact GEM correspondence

For Human-GEM and Mouse-GEM, the package no longer identifies a nutrient by
searching arbitrary reaction text alone. It now:

1. identifies exchange reactions from reaction roles;
2. reads each exchange column from the stoichiometric matrix;
3. resolves its unique extracellular metabolite;
4. joins the metabolite to `gem$metabolite_meta`;
5. compares normalized **exact aliases** with the GEM metabolite name;
6. records the matched GEM metabolite ID, name, and match method.

A reaction description that merely contains a nutrient name is therefore not
enough to open the exchange. Pattern matching remains only as a compatibility
fallback for generic user-supplied GEMs that lack metabolite metadata.

If two preset compounds resolve to the same official GEM exchange, construction
fails instead of silently choosing one.

## Deliberate exclusions

Lipoprotein-bound cholesterol, phospholipids, triacylglycerols, and
albumin-bound fatty acids are not opened by this change. Total plasma abundance
does not establish a freely available extracellular uptake bound, and different
cell types access these pools through distinct transport, lipolysis, and
receptor-mediated mechanisms. Such nutrients should be supplied through a
measured or explicitly assumed custom medium.

## References

- Cantor JR et al. Physiologic Medium Rewires Cellular Metabolism and Reveals
  Uric Acid as an Endogenous Inhibitor of UMP Synthase. *Cell* (2017).
  DOI: 10.1016/j.cell.2017.03.023.
- Psychogios N et al. The Human Serum Metabolome. *PLoS ONE* (2011).
  DOI: 10.1371/journal.pone.0016957.
- Gardner GL and Stuart JA. Tumor microenvironment-like conditions alter
  pancreatic cancer cell metabolism and behavior. *Am J Physiol Cell Physiol*
  (2024). DOI: 10.1152/ajpcell.00452.2024.
- Sullivan MR et al. Quantification of microenvironmental metabolites in murine
  cancers reveals determinants of tumor nutrient availability. *eLife* (2019).
  DOI: 10.7554/eLife.44235.
