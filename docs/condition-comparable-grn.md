# Pando condition-aware GRN contract

This page describes the implemented Pando–RegCompass interface. Equations are
centralized in [Mathematical model](mathematical-model.md).

## Scope

`rc_regcompass_step_grn()` fits `Pando::infer_condition_grn()` separately for
each broad cell type while using all eligible conditions within that type.
Biological sample metadata are not model inputs or fitting gates.

The canonical configuration uses:

```r
pando_infer_args = list(
  candidate_screen = "motif_domain",
  condition_mix = 0.5,
  condition_weight = "equal",
  outer_nfolds = 5L,
  inner_nfolds = 5L,
  scale = TRUE
)
```

## Implemented model behavior

- candidate edges are defined by motif and regulatory-domain support;
- the fitted predictor combines TF RNA and peak accessibility;
- conditions are fitted jointly with a sparse-group multitask objective;
- coefficients may differ in magnitude, sign, or active status by condition;
- transforms and lambda selection are learned inside each outer training fold;
- a support-constrained common-metric refit is retained for interpretation;
- full-data projections are not used for the primary penalty.

The condition columns are not independent elastic-net fits.

## Candidate screening

RegCompass requires `candidate_screen = "motif_domain"`. The wrapper sets this
value when omitted and stops if another value is supplied.

Pando may expose additional candidate-screening modes, but they are not accepted
by the current RegCompass Stage 1 interface.

## `ConditionGRNFit` v5 fields

Important fields include:

```text
edge_table
beta_condition
beta_shared
delta_condition
estimability_mask
support_mask
active_mask
predictor_transform
response_transform
target_fit
condition_rsq_oof
target_rsq_oof_pooled
projection_common_oof
projection_condition_full_oof
projection_global_common_oof
cell_provenance
```

Unavailable fit coefficients remain `NA`; an estimable inactive coefficient
remains numeric zero. The public and persisted contract contains no stored
baseline-condition coefficient, contrast matrix, or comparison mask.

## Absolute condition effects

`beta_condition` is the absolute condition effect on the shared
within-cell-type equal-condition coordinate. RegCompass uses this view for
supported-gene selection and the primary Stage 4 projection.

Conditions may have different active edges and opposite coefficient directions.
They remain comparable because the candidate dictionary and coefficient scale
are shared inside the fitted broad cell type.

## Primary projection contract

The canonical Stage 4 path requires:

```text
projection_component = "condition"
origin = "oof"
support_policy = "pairwise_common" or "global_common"
```

- two-condition analyses use pairwise-common support;
- analyses with more than two conditions use global-common support;
- Pando computes the single-cell score before SuperCell aggregation;
- a non-estimable edge contributes exactly zero at the projection-contribution
  layer;
- the structural zero enters target summation, metacell averaging, GPR
  aggregation, reaction expression and the main penalty;
- structural-zero masks and fractions remain available for audit;
- RegCompass averages completed target-gene scores by exact membership;
- pooled OOF reliability and one pooled broad-cell-type calibration scale are
  used.

When a target-level regulatory modifier remains unavailable or non-finite,
RegCompass uses a neutral modifier and therefore returns exactly the RNA-only
support for that gene–metacell entry. The fallback is explicitly annotated.

## Stage ownership

RegCompass forwards:

```text
pando_initiate_args → Pando::initiate_grn()
pando_motif_args    → Pando::find_motifs()
pando_infer_args    → Pando::infer_condition_grn()
```

RegCompass controls the object, assay names, motif object, genome, metadata
columns, GEM target genes, network name, minimum condition size, error policy,
and `BPPARAM`. Nested overrides of these fields are rejected.

`aggregate_rna_col` and `aggregate_peaks_col` are also rejected because Stage 1
fits paired single cells and Stage 2 owns metacell construction.

## Persisted artifacts

Stage 1 writes absolute-condition artifacts only:

```text
pando_group_status.tsv.gz
pando_tf_peak_gene_condition_all.tsv.gz
pando_tf_peak_gene_condition_active.tsv.gz
pando_tf_peak_gene_condition_effect_all.tsv.gz
pando_tf_peak_gene_condition_effect_active.tsv.gz
pando_tf_peak_gene_universal.tsv.gz
pando_condition_network_index.tsv.gz
pando_condition_fit_diagnostics.tsv.gz
pando_edge_predictor_transforms.tsv.gz
pando_condition_grn_fits.rds
pando_objects/condition_grn_fit_v5.rds
```

The legacy `condition_effect` filenames are retained for file compatibility, but
their numeric definition is the absolute condition coefficient and the tables
contain `effect_definition = "absolute_condition_coefficient"`.

## Genome-build requirement

Human analyses may use the bundled hg38 regulatory-region union. Mouse analyses
must supply a species- and build-matched `GRanges` object through
`pando_initiate_args$regions`. The region build must match both the ATAC peak
coordinates and the genome used for motif scanning.
