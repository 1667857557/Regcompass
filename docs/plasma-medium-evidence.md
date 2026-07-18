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

`mouse_plasma` uses the murine plasma and tumor-interstitial-fluid metabolite
detection evidence reported by Sullivan et al. Concentrations that were not
validated as directly comparable between human and mouse are retained as
availability-only rows (`concentration_mM = NA`) rather than borrowing the human
HPLM value.

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

The concentration column is provenance. Except for explicitly flagged
glucose/lactate/glutamine sensitivity scenarios, RegCompass does not interpret a
medium concentration as a measured uptake flux.

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
- Sullivan MR et al. Quantification of microenvironmental metabolites in murine
  cancers reveals determinants of tumor nutrient availability. *eLife* (2019).
  DOI: 10.7554/eLife.44235.
