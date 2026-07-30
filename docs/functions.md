# RegCompassR public API

See [Workflow](workflow.md), [run modes](run-modes-and-stepwise-workflow.md),
[stage contracts](stage-interface-contracts.md), and the
[metacell graph contract](metacell-graph-contract.md).

## Complete workflows

| Function | Purpose |
|---|---|
| `rc_run_regcompass()` | Run all six stages with automatic standard/condition-aware Pando routing. |
| `rc_run_regcompass_one_shot()` | Convenience wrapper around the complete workflow. |

`condition_col` may be `NULL`, absent from metadata, single-level, or
multi-level. The selected route is returned in `result$analysis_mode`.

## Restartable stages

| Stage | Function | Main output |
|---:|---|---|
| 1 | `rc_regcompass_step_grn()` | standard Pando networks or canonical condition-aware fits |
| 2 | `rc_regcompass_step_metacells()` | cell-type-independent graphs, condition-pure membership, and aggregated RNA/ATAC counts |
| 3 | `rc_regcompass_step_meta_modules()` | supported genes, complete-GPR cores, reaction catalogue |
| 4 | `rc_regcompass_step_layer1()` | RNA support, regulatory projection, GPR reaction expression |
| 5 | `rc_regcompass_step_layer2()` | shared medium-specific models and directional LP scores |
| 6 | `rc_regcompass_step_results()` | annotations, rankings, summaries, and available contrasts |

## Stage 1 routing

```r
step1 <- rc_regcompass_step_grn(
  object,
  gem,
  outdir,
  genome,
  condition_col = "condition",
  celltype_col = "cell_type",
  pando_args = list(
    min_cells = 100L,
    min_model_rsq = 0.1,
    pando_infer_args = list(...)
  )
)
```

Nested arguments are routed to:

| Bundle | Target |
|---|---|
| `pando_initiate_args` | `Pando::initiate_grn()` |
| `pando_motif_args` | `Pando::find_motifs()` |
| `pando_infer_args` | selected standard or condition-aware Pando inference function |

In condition mode, `candidate_screen="motif_domain"`,
`condition_weight="equal"`, and `scale=TRUE` are required. In standard mode,
original `Pando::infer_grn()` is used and no condition coefficient is produced.

Important outputs:

```r
step1$params$analysis_mode
step1$grn_result$condition_coefficients_calculated
step1$grn_result$standard_pando_objects   # standard mode
step1$grn_result$condition_grn_fits       # condition mode
step1$grn_result$tf_peak_gene_condition
```

## Stage 2 SuperCell graph arguments

```r
metacell_args = list(
  rna_reduction = "pca",
  rna_dims = 1:30,
  atac_reduction = "lsi",
  atac_dims = 2:30,
  gamma = 30L,
  k.knn = 5L,
  seed = 12345L
)
```

RegCompass constructs a cell-type-specific multimodal embedding and calls
`SuperCell::SCimplify_by_graph_group_from_embedding()` with:

```text
cell.graph.group = broad cell type
cell.split.condition = condition or NULL
```

`cell.graph.group` creates one independent kNN graph per cell type.
`cell.split.condition` acts only after graph clustering, so all conditions share
the cell-type graph but final metacells remain condition-pure. Embedding blocks
are standardized within cell type across all conditions and divided by the
square root of their dimension count before concatenation.

It does not use `sample`, expose a combined stratum column, or build separate
condition graphs. Cache provenance is stored in
`step2$pooled$cache_contract`; the complete formal contract is stored in
`step2$pooled$input_design`.

## Stage 4 regulatory support

```r
step4 <- rc_regcompass_step_layer1(
  grn,
  metacells,
  meta_modules,
  gem,
  outdir,
  comparison_support = "auto",
  regulatory_alpha = 1,
  gpr_and_method = "min"
)
```

- condition mode uses common-support outer-heldout projections;
- standard mode uses standard Pando coefficients and records full-fit origin;
- both calculate TF RNA × peak ATAC per cell before metacell aggregation;
- a non-finite target modifier uses neutral `R=0` and equals RNA-only support;
- `regulatory_alpha` is fixed at `1`;
- GPR AND accepts `min`, `median`, or `mean`; OR branches are additive;
- missing final reaction expression is converted to `E=0` and penalty `1`.

## GEM and medium

| Function | Purpose |
|---|---|
| `rc_prepare_gem()` | Load and validate a pinned Human-GEM or Mouse-GEM. |
| `rc_make_medium_scenarios()` | Build preset or custom exchange bounds. |
| `rc_build_reaction_annotations()` | Build reaction names, formulas, GPRs, and cross-references. |
| `rc_attach_reaction_annotations()` | Attach annotations to an existing result. |

## Condition analysis

| Function | Purpose |
|---|---|
| `rc_test_condition_reactions()` | Compare fixed reaction-direction targets when multiple conditions exist. |
| `rc_report_condition_directions()` | Summarize forward and reverse targets. |
| `rc_plot_condition_reaction()` | Plot one reaction direction across conditions. |
| `rc_select_gene_reactions()` | Select scored reactions by metabolic gene. |

For one condition, use `result$reaction_ranking` and
`result$condition_summary`; `result$condition_contrast` is empty.

Metacell tests are within-dataset inference, not biological-replicate inference.
