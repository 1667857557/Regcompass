# Condition-comparability invariants

Equations are in [Tutorial 3](tutorial-03-mathematical-model.md). This page lists
only structural validation rules.

## Pando fit contract

RegCompass requires the canonical unversioned `pando_condition_grn_fit` with:

- one fitted broad cell type;
- a shared candidate supergraph and equal-condition coordinate;
- nested outer-heldout projections with exactly-once cell assignment;
- `coefficient_estimable_mask`;
- `projectable_structural_zero_mask`;
- `projection_support_mask`;
- finite condition-full and common-support OOF target scores;
- no full-data condition projection in the primary penalty.

Unavailable coefficients remain `NA`; estimable inactive coefficients remain
numeric zero. A non-estimable edge contribution is fixed at zero in that
condition.

## Primary projection

`condition_full_oof` is primary. Jointly estimable edges form the common-support
component. Their difference is the condition-unique component.

An exact-zero predictor remains in the shared candidate supergraph, receives no
fitted coefficient, and contributes zero. The zero decision for held-out cells
uses only the corresponding training fold.

## Candidate and stage ownership

Stage 1 requires `candidate_screen = "motif_domain"`,
`condition_weight = "equal"`, and `scale = TRUE`.

```text
pando_initiate_args → Pando::initiate_grn()
pando_motif_args    → Pando::find_motifs()
pando_infer_args    → Pando::infer_condition_grn()
```

Stage 2 owns metacell construction; Pando aggregation columns are rejected.

## Metacell invariants

- one independent graph per broad cell type;
- all conditions jointly participate in that cell-type graph;
- condition is applied after graph clustering;
- each input cell maps to exactly one metacell;
- final metacells contain one cell type and one condition;
- RNA and ATAC metacell matrices have identical ordered IDs;
- no sample-derived or combined condition-by-cell-type grouping is created.

## Shared metabolic-model invariants

Within one medium, all compared units share reaction order, stoichiometry,
bounds, target direction, target-flux fraction and target-specific `vmax`. A
mismatch is an error rather than a condition difference.

The retired depth-matching, common-depth, alpha-sensitivity, zero-support-
sensitivity and link-saturation-propagation branches are not stage invariants and
are not persisted.

Public API: [functions.md](functions.md).
