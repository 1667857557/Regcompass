# Canonical defaults

The current main workflow uses:

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
  rna_reduction = "pca",
  atac_reduction = "lsi"
)

layer1_args = list(
  gpr_and_method = "min"
)
```

Each broad cell type receives one independent multimodal graph. All conditions
of that cell type participate jointly in the graph, and condition-pure metacells
are assigned after graph clustering.

`condition_full_oof` is the primary Layer 1 and Layer 2 route. The
component retained for decomposition. It does not replace the primary route.

A non-estimable edge side has an unavailable coefficient and a projectable
`R = 0`, exactly recovering RNA-only support.

The canonical schema does not calculate depth matching, common-depth
restriction, alpha sensitivity, zero-support sensitivity, or link-saturation
propagation.

Equations: [Tutorial 3](tutorial-03-mathematical-model.md). Public API:
[functions.md](functions.md).
