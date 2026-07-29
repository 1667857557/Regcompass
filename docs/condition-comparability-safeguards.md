# Condition-comparability safeguards

Mathematical definitions are in [Mathematical model](mathematical-model.md).
This page lists validation rules that affect analysis.

## Pando fit requirements

RegCompass requires a Pando `ConditionGRNFit v5` with:

- one fitted broad cell type;
- complete condition labels and cell provenance;
- equal-condition coefficient scaling;
- nested outer-heldout cell projections;
- explicit estimability, support, activity, and comparison masks;
- complete OOF assignment and coverage;
- no full-data projection in the primary penalty path.

Unavailable coefficients must be `NA`; estimable inactive coefficients remain
zero.

## Primary comparison support

- two-condition analyses use pairwise-common estimability;
- multi-condition analyses use global-common estimability;
- condition-estimable and strict projections are diagnostic only;
- reference-condition effects are interpretation outputs;
- condition-specific edges do not enter the primary common-support penalty.

## Candidate screening

`candidate_screen = "motif_domain"` is the canonical mode.

`pooled_within_condition` is a response-dependent sensitivity mode and is
ineligible for primary penalty construction.

## Stage ownership

RegCompass routes:

```text
pando_initiate_args → initiate_grn
pando_motif_args    → find_motifs
pando_infer_args    → infer_condition_grn
```

It rejects nested overrides of managed objects, assays, genome, metadata
columns, target genes, network name, minimum condition size, error policy, and
`BPPARAM`.

Pando aggregation columns are rejected because Stage 2 owns metacell
construction.

## Parallel execution

| `parallel` | `BPPARAM` | Stage 1 route |
|---|---|---|
| `FALSE` | any | serial |
| `TRUE` | `BiocParallelParam` | supplied backend |
| `TRUE` | `NULL` or `FALSE` | Pando native map |
| any | `TRUE` | error |

The resolved route is stored in `step1$params$pando_parallel`.

## Metacell checks

- strata are condition × broad cell type;
- every input cell maps to exactly one metacell;
- RNA and ATAC metacell matrices must contain identical ordered IDs;
- mixed condition or mixed broad-cell-type metacells are rejected;
- cache reuse requires identical cells, assays, reductions, dimensions, seed,
  gamma, and thresholds.

## Shared metabolic-model checks

Within one medium comparison, all units must share:

- reaction order;
- stoichiometric matrix and checksum;
- bounds;
- target direction;
- target-flux fraction;
- target-specific `vmax`.

A mismatch is an error rather than a reported condition difference.

## Genome build

Human analyses may use the bundled hg38 regulatory regions. Mouse analyses must
supply a build-matched region set. The region build must agree with both ATAC
coordinates and the motif-scanning genome.
