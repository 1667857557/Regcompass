# Condition-comparability safeguards

## Why this change is required

Pando estimates condition coefficients on one shared TF–peak–target dictionary. An edge can nevertheless be non-estimable in one condition because the TF, peak, or final TF-by-peak predictor has no eligible variation in that condition.

A coefficient fixed to zero by the eligibility mask is therefore not equivalent to an estimated zero coefficient. For an explicit reference condition, RegCompass now requires

```text
comparable_to_reference[e, c] =
    eligible[e, c] && eligible[e, reference]
```

before an edge-level condition effect can enter the regulatory modifier or penalty calculation.

The complete edge tables retain all rows and add `comparable_to_reference`. The active condition-effect table excludes non-comparable rows.

## Candidate-edge policy

The fitted Pando predictor is

```text
TF_RNA * peak_ATAC
```

Marginal TF–target and peak–target correlations can both be near zero when the interaction predictor is strongly associated with the target. RegCompass therefore defaults to

```r
candidate_screen = "motif_domain"
```

which retains the structurally supported motif/domain edge dictionary and lets the shared elastic-net model perform coefficient selection. Users can still explicitly request `condition_union` or `pooled`, but those faster modes impose marginal-correlation assumptions that are not guaranteed by the interaction model.

## Genome-build safety

The Pando fork bundles an hg38 conserved-element region set. It is invalid for mouse ATAC coordinates. Mouse analyses must supply a species- and build-matched `GRanges` object through

```r
pando_initiate_args = list(regions = mouse_regions)
```

The region genome build must match both the ATAC peak coordinates and the genome object passed to motif scanning. RegCompass now stops rather than silently applying hg38 regions to mouse input.

## Unchanged model properties

This correction does not change:

- the shared TF–peak–target edge coordinate system;
- pooled predictor and target transformations;
- target-specific shared lambda paths;
- independently estimated condition coefficient columns at a fixed lambda;
- the shared GEM and stoichiometric reaction space.
