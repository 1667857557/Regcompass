# Condition-comparability safeguards

## Explicit comparison support

Pando fits condition coefficients on one shared TF–peak–target dictionary. An
edge can still be non-estimable in one condition because its TF, peak, or final
`TF RNA × peak ATAC` predictor lacks eligible variation.

A coefficient fixed to zero by an eligibility constraint is not equivalent to
an eligible coefficient estimated as zero. Pando 1.5.0 therefore exports:

```text
comparison_mask[e, c] =
  eligibility_mask[e, c] && eligibility_mask[e, reference]
```

RegCompass validates and requires this explicit matrix. It no longer silently
reconstructs comparison support for older Pando objects. The complete condition
effect table retains every edge and adds `comparable_to_reference`; the active
effect table excludes non-comparable rows before Layer 1 or penalty calculation.

## Candidate-edge policy

The fitted Pando predictor is:

```text
TF_RNA × peak_ATAC
```

Marginal TF-target and peak-target correlations can both be near zero even when
the interaction predictor is informative. RegCompass therefore defaults to:

```r
candidate_screen = "motif_domain"
```

This retains the structurally supported motif/domain dictionary and lets the
shared elastic-net model select coefficients. Users can explicitly request
`pooled_within_condition` as a marginal-screen sensitivity analysis. Its response-dependent screen makes its projection ineligible for penalty construction.

## Pando argument ownership

RegCompass routes:

```text
pando_initiate_args → initiate_grn
pando_motif_args    → find_motifs
pando_infer_args    → infer_condition_grn
```

It rejects nested overrides of the object, assays, motif object, genome,
condition/cell-type columns, GEM target genes, network name, minimum condition
size, error policy, and `BPPARAM`. This prevents a nested argument from changing
the fit coordinate system while downstream extraction still assumes the
RegCompass-managed contract.

`aggregate_rna_col` and `aggregate_peaks_col` are rejected in canonical Stage 1.
Pando is fitted on paired single cells; RegCompass Stage 2 owns metacell
aggregation.

## Parallel routing

For `rc_regcompass_step_grn()`:

```text
parallel = FALSE                         → serial
parallel = TRUE + BiocParallelParam      → supplied backend
parallel = TRUE + BPPARAM NULL/FALSE     → Pando native map
BPPARAM = TRUE                           → error
```

The resolved route is stored in `step1$params$pando_parallel`. The one-shot
workflow uses `upstream_workers` and a stage-scoped backend, so nested
`pando_infer_args$parallel` should not be supplied.

## Genome-build safety

Pando bundles hg38 conserved-element/SCREEN regions. Those coordinates are
invalid for mouse ATAC peaks. Mouse analyses must supply a species- and
build-matched `GRanges` object through:

```r
pando_initiate_args = list(regions = mouse_regions)
```

The region genome build must match both the ATAC peak coordinates and the genome
object passed to motif scanning. RegCompass stops rather than silently applying
hg38 regions to mouse input.

## Unchanged model properties

These safeguards do not change:

- the shared TF–peak–target coordinate system;
- pooled predictor and target transformations;
- target-specific shared lambda paths;
- condition-sparse selection followed by common-metric refit;
- the shared GEM and stoichiometric reaction space.
