# Canonical defaults

The canonical RegCompass settings used by the current workflow are:

```r
pando_infer_args = list(
  candidate_screen = "motif_domain",
  peak_cor = 0,
  condition_mix = 0.5,
  condition_weight = "equal",
  scale = TRUE
)

metacell_args = list(
  gamma = 30L,
  depth_balance = FALSE
)

layer1_args = list(
  projection_component = "condition",
  comparison_support = "auto",
  regulatory_alpha = 1,
  gpr_and_method = "min"
)
```

One fixed gamma is used across every condition × broad-cell-type stratum.
Non-estimable regulatory edge contributions are structural zeros in the main
analysis. A non-finite target-level modifier falls back to RNA-only support.
