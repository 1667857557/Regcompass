# RegCompassR workflow

This page describes stage responsibilities, inputs, and outputs. Mathematical
definitions are centralized in [Mathematical model](mathematical-model.md).

## Data flow

```text
paired RNA+ATAC cells
→ condition-aware Pando GRNs
→ condition × broad-cell-type metacells
→ supported metabolic genes and core reactions
→ merged reaction catalogue
→ RNA support modified by regulatory scores
→ GPR-based reaction penalties
→ shared medium-specific metabolic model
→ directional LP scores
→ condition comparisons
```

## Stage 1: `rc_regcompass_step_grn()`

**Purpose:** fit one condition-aware Pando model per broad cell type.

**Main inputs:** paired RNA+ATAC Seurat object, GEM, genome, condition column,
cell-type column, Pando arguments.

**Key behavior:**

- targets are GEM GPR genes present in the RNA assay;
- `candidate_screen = "motif_domain"` is the default;
- outer-heldout projections are generated within each fitted cell type;
- absolute condition coefficients are used downstream;
- reference-condition effects are interpretation outputs;
- mouse analyses require build-matched regulatory regions.

**Main outputs:**

```r
step1$grn_result$condition_grn_fits
step1$grn_result$condition_fit_status
step1$grn_result$tf_peak_gene_condition
step1$grn_result$tf_peak_gene_condition_effect
```

Detailed model fields and equations: [Pando condition contract](condition-comparable-grn.md)
and [Mathematical model](mathematical-model.md).

## Stage 2: `rc_regcompass_step_metacells()`

**Purpose:** construct multimodal SuperCell metacells.

**Grouping:** condition × broad cell type. Metacells mixing either field are
rejected.

**Important parameters:** RNA/ATAC reductions and dimensions, `gamma`, `seed`,
`min_cells_per_stratum`, `min_metacell_size`, and
`min_metacells_per_stratum`.

**Main outputs:** metacell RNA and ATAC counts, membership table, metacell
metadata, cache contract, and stratum status.

Changing cells, assays, reductions, dimensions, seed, gamma, or thresholds
invalidates the cache.

## Stage 3: `rc_regcompass_step_meta_modules()`

**Purpose:** map active condition coefficients to metabolic genes and reaction
sets.

**Processing:**

1. identify supported metabolic genes;
2. require a complete GPR branch for a core reaction;
3. add reactions from the core subsystem;
4. add direct KEGG/Reactome equivalents;
5. add direct master-Rhea equivalents;
6. merge condition- and cell-type-specific sets into one reaction catalogue.

Stage 3 does not run FASTCORE and does not create the final GEM.

**Main outputs:** supported genes, core gene-reaction mapping, reaction
membership, merged core reactions, and merged reaction catalogue.

## Stage 4: `rc_regcompass_step_layer1()`

**Purpose:** calculate metacell gene support and reaction penalties.

**Key behavior:**

- RNA support is estimated from metacell counts;
- the primary regulatory input is Pando's outer-heldout common-support score;
- single-cell regulatory scores are averaged using exact SuperCell membership;
- reference-condition coefficient effects are not used for the primary penalty;
- `gpr_and_method = "min"` is the default;
- `"median"` and `"mean"` are available for sensitivity analysis;
- RNA-only, condition-full, depth, and alpha routes are diagnostic or sensitivity outputs.

**Main outputs:** gene support, regulatory modifiers, reaction support, reaction
penalties, and route diagnostics.

## Stage 5: `rc_regcompass_step_layer2()`

**Purpose:** build a shared medium-specific model and score reaction directions.

**Key behavior:**

- one global FASTCORE completion is run per medium;
- the completed model is reused for every condition and metacell in that medium;
- forward and reverse directions are scored separately;
- lower normalized penalty indicates stronger network-constrained support;
- the score is not a measured flux.

**Important parameters:** medium scenario, `target_direction`, `omega`, solver,
`completion_time_limit`, `fastcore_epsilon`, `max_support_reactions`, and
`strict`.

**Main outputs:** model cache, model diagnostics, directional penalties, `vmax`,
feasibility, and comparison tables.

Medium choices: [Medium presets](medium-presets.md).

## Stage 6: `rc_regcompass_step_results()`

**Purpose:** assemble annotations, rankings, evidence, provenance, and condition
contrasts.

Condition comparisons must use the same reaction, direction, medium, broad cell
type, model, bounds, and target-flux fraction.

Metacell P values describe within-dataset condition-associated separation. They
are not sample- or donor-level biological-replicate inference.

## Optional targeted scoring

`rc_regcompass_step_target_union()` scores direct KEGG, Reactome, or master-Rhea
equivalents of selected reaction anchors in the existing Stage 5 model. It does
not rebuild the model or rerun FASTCORE.

## Restart boundaries

| Earliest stage to rerun | Changes |
|---|---|
| Stage 1 | RNA/ATAC data, cell labels, genome, regions, motifs, Pando fitting or filtering arguments |
| Stage 2 | metacell reductions, dimensions, gamma, seed, thresholds, or membership inputs |
| Stage 3 | GPR rules, subsystem table, or reaction cross-references |
| Stage 4 | regulatory alpha, RNA support settings, GPR aggregation rule, or metacell evidence |
| Stage 5 | medium, bounds, target direction, omega, solver, or model-completion settings |
| Stage 6 | annotations, reporting filters, or contrast settings only |

Public API index: [functions.md](functions.md).
