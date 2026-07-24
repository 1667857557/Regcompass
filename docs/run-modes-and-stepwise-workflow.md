# RegCompassR 1.8.2 tutorial index

Choose the lowest level that provides the required control. All levels use the same GRN-first biological model and defaults.

| Level | Use | Tutorial |
|---|---|---|
| 1 | validated one-shot analysis | [Quick start](tutorial-01-quick-start.md) |
| 2 | explicit stages and continuation gates | [Stepwise audit](tutorial-02-stepwise-audit.md) |
| 3 | restart, sensitivity, resource allocation, and failure diagnosis | [Advanced restart](tutorial-03-advanced-restart.md) |

The exact input/output requirements for every stage are summarized in [Stage input-output contracts](stage-interface-contracts.md).

## Canonical workflow

```text
single-cell RNA normalization
→ cell-type-shared ATAC TF-IDF across conditions
→ Pando GRN per condition × cell type (peak_cor = 0.01)
→ cell-type-guided SuperCell2 metacells within condition (gamma = 75)
→ complete-GPR core reactions
→ subsystem + KEGG/Reactome + master-Rhea expansion
→ local FASTCORE feasibility completion
→ global union GEM
→ integrated RNA+ATAC reaction expression
→ directional COMPASS-like LP scoring
→ reaction annotation and condition comparisons
```

Stage 3-6 validate stage classes, workflow metadata, GEM fingerprints, core-target provenance, and ordered scoring units before connecting objects.

## Optional analyses after Layer 2

- [Expanded target scoring](target-union-scoring.md): select previous core reactions or GPR genes, expand their annotation-related context, and score every expanded reaction in the exact cached union GEM from the original run.
- [Condition-associated reaction statistics](condition-reaction-statistics.md): compare the same reaction, direction, medium, and cell type across conditions.

## Parallel units

| Stage | Parallel unit |
|---|---|
| GRN | condition × cell type |
| Metacells | no workflow-level BiocParallel loop |
| Meta-modules | local FASTCORE completion per module |
| Layer 1 | GPR/reaction capacity |
| Layer 2 | shared model × metacell |
| Expanded targets | reused union model × metacell |
| Results | serial assembly |

On Linux use `BiocParallel::MulticoreParam` or `parallel_backend = "multicore"`. Keep Pando's inner `parallel = FALSE` and set numerical-library thread counts to one.

## Required input

- paired-cell RNA and ATAC counts in a Seurat object;
- a Signac `ChromatinAssay` and matching genome coordinates;
- complete condition and cell-type metadata;
- `Pando::motifs` or another compatible PFM/PWM collection;
- a validated human or mouse GEM;
- an installed LP solver.

Sample metadata are optional provenance and are not used for balancing or grouping. Metacell-level condition comparisons are descriptive within-dataset analyses unless independent biological replication is supplied.
