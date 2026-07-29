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
  reference_condition = "Control",
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

`candidate_screen = "motif_domain"` is the canonical mode. It avoids
response-dependent marginal screening before nested validation.

`candidate_screen = "pooled_within_condition"` remains available as a
sensitivity mode. Its projection is marked ineligible for primary penalty
construction because candidate screening uses the response outside the nested
OOF loop.

## `ConditionGRNFit` v5 fields

Important fields include:

```text
edge_table
beta_condition
beta_shared
delta_condition
contrast
estimability_mask
support_mask
active_mask
comparison_mask
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

Unavailable coefficients remain `NA`; an estimable inactive coefficient remains
numeric zero.

## Absolute coefficients and reference effects

RegCompass retains two coefficient views:

- **absolute condition coefficients:** used to identify supported metabolic genes and to generate the primary Stage 4 projection;
- **reference-condition effects:** used in interpretation tables only.

The effect table includes `comparable_to_reference`, derived from the explicit
Pando comparison mask. RegCompass validates the mask and does not reconstruct it
from older objects.

## Primary projection contract

The canonical Stage 4 path requires:

```text
projection_component = "condition"
origin = "oof"
support_policy = "pairwise_common" or "global_common"
```

- two-condition analyses use pairwise-common estimability;
- analyses with more than two conditions use global-common estimability;
- condition-estimable and strict projections are diagnostic only;
- the primary projection uses absolute condition coefficients, not reference contrasts;
- Pando computes the single-cell score before SuperCell aggregation;
- RegCompass averages the completed target-gene scores by exact membership;
- pooled OOF reliability and one pooled broad-cell-type calibration scale are used.

Condition-specific edges remain available in network and diagnostic outputs but
do not enter the primary common-support penalty.

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

Stage 1 writes:

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

## Genome-build requirement

Human analyses may use the bundled hg38 regulatory-region union. Mouse analyses
must supply a species- and build-matched `GRanges` object through
`pando_initiate_args$regions`. The region build must match both the ATAC peak
coordinates and the genome used for motif scanning.
