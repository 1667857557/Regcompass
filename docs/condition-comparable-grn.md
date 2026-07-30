# Pando condition-aware GRN contract

This page describes the implemented Pando–RegCompass interface. Equations are
centralized in [Tutorial 3](tutorial-03-mathematical-model.md).

## Scope

`rc_regcompass_step_grn()` calls `Pando::infer_condition_grn()` separately for
each broad cell type while using all eligible conditions of that type. Sample
metadata are not model inputs or fitting gates.

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

Conditions share the candidate supergraph, coefficient units and fold-local
transforms. Coefficients may differ in magnitude, sign and active status.

## Canonical fit fields

```text
edge_table
beta_condition_std
coefficient_estimable_mask
projectable_structural_zero_mask
projection_support_mask
support_mask
active_mask
predictor_transform
response_transform
condition_rsq_oof
target_rsq_oof_pooled
projection_condition_full_oof
projection_common_oof
projection_global_common_oof
```

An unavailable coefficient remains `NA`; an estimable inactive coefficient is
numeric zero. A projectable structural zero has no fitted coefficient and a
fixed projection contribution of zero.

## Primary projection contract

The canonical Stage 4 call is:

```r
Pando::project_condition_grn_primary_cells(
  object = pando_grn_data,
  fit = fit,
  scale = "std"
)
```

The primary score is condition-full outer-heldout projection. For each condition:

- estimable edges contribute their outer-fold condition coefficients;
- jointly estimable edges form the common-support component;
- unilateral edges contribute only in their estimable condition;
- non-estimable sides contribute zero;
- bilaterally non-estimable edges contribute zero in both conditions;
- exact-zero predictors remain represented in the candidate supergraph.

The condition-unique component is the condition-full projection minus the
common-support projection. Projection occurs on paired single cells before exact
SuperCell aggregation. No coefficient, centre or scale is refitted after
aggregation.

## RegCompass fields

Stage 1 adds auditable edge-level columns:

```text
coefficient_estimable
projectable_structural_zero
projection_supported
```

Stage 4 stores:

```text
gene_projection_condition_full_oof
gene_projection_common_oof
gene_projection_condition_unique_oof
reaction_expression_condition_full_oof
reaction_expression_common_oof
```

Stage 5 stores:

```text
penalty_condition_full_oof
penalty_common_oof
penalty_condition_unique_increment
penalty_rna_only
```

All Stage 5 routes use the same GEM, bounds, reaction order, direction and
`vmax`.

## Stage ownership

```text
pando_initiate_args → Pando::initiate_grn()
pando_motif_args    → Pando::find_motifs()
pando_infer_args    → Pando::infer_condition_grn()
```

RegCompass controls assay names, metadata columns, target genes, network name,
minimum condition size, error policy and parallel backend. Stage 2 owns
metacell construction.

## Persisted Stage 1 artifacts

```text
pando_group_status.tsv.gz
pando_tf_peak_gene_condition_all.tsv.gz
pando_tf_peak_gene_condition_active.tsv.gz
pando_tf_peak_gene_universal.tsv.gz
pando_condition_network_index.tsv.gz
pando_condition_fit_diagnostics.tsv.gz
pando_edge_predictor_transforms.tsv.gz
pando_condition_grn_fits.rds
```

Human analyses may use bundled hg38 regions. Mouse analyses must provide a
species- and build-matched `GRanges` object.

Public API: [functions.md](functions.md).
