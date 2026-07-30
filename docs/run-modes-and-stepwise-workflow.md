# Run modes and stepwise workflow

RegCompass uses one public workflow and resolves the Pando mode automatically.

## Mode selection

```text
condition has >= 2 levels  -> condition_grn
condition has 1 level      -> standard_pando
condition omitted/absent   -> standard_pando
```

### Condition mode

`Pando::infer_condition_grn()` fits one canonical `pando_condition_grn_fit`
contract per broad cell type. The primary penalty uses condition-full
outer-heldout TF RNA × peak ATAC projection. Jointly estimable edges form the
common-support component; every non-estimable edge side is a projectable
structural zero.

### Standard mode

`Pando::infer_grn()` runs within each broad cell type. Standard coefficients are
projected to paired cells and aggregated to metacells. No condition coefficients
or condition contrast are calculated.

## One-shot workflow

```r
result <- rc_run_regcompass(
  object = object,
  gem = gem,
  outdir = "result",
  genome = genome,
  condition_col = "condition",
  celltype_col = "cell_type"
)
```

`result$analysis_mode` reports the selected route.

## Stepwise workflow

```r
step1 <- rc_regcompass_step_grn(...)
step2 <- rc_regcompass_step_metacells(...)
step3 <- rc_regcompass_step_meta_modules(...)
step4 <- rc_regcompass_step_layer1(...)
step5 <- rc_regcompass_step_layer2(...)
result <- rc_regcompass_step_results(...)
```

Stage 1 and Stage 2 independently resolve the design and their `analysis_mode`
values must agree.

## SuperCell graph and purity inputs

Stage 2 calls `SCimplify_by_graph_group_from_embedding()` with:

```text
cell.graph.group       = broad cell type
cell.split.condition   = condition, or NULL when omitted
gamma                   = requested graining level
```

`cell.graph.group` partitions cells before neighbour construction. All conditions
within a cell type share one standardized multimodal embedding and graph.
`cell.split.condition` is applied after graph clustering, so final metacells are
condition-pure without separate condition graphs.

No combined metadata stratum or sample-derived grouping is created. Provenance
records `one_independent_graph_per_cell_type`,
`all_conditions_joint_within_cell_type_graph`, and
`temporary_combined_stratum = FALSE`.

## Restart boundaries

- motifs, regions, targets or Pando fitting changed: rerun Stage 1 onward;
- reductions, dimensions, gamma or SuperCell settings changed: rerun Stage 2 onward;
- GPR or catalogue annotations changed: rerun Stage 3 onward;
- projection, RNA support or GPR aggregation changed: rerun Stage 4 onward;
- medium, union GEM, FASTCORE or LP controls changed: rerun Stage 5 onward.

## Result behavior

Multiple conditions produce primary condition-full summaries and pairwise
contrasts. Common-support and RNA-only routes remain decomposition/control
outputs. A single effective condition produces rankings and summaries with an
empty condition contrast.

See [Tutorial 1](tutorial-01-quick-start.md),
[Tutorial 2](tutorial-02-stepwise-audit.md),
[Tutorial 3](tutorial-03-mathematical-model.md),
[Tutorial 4](tutorial-04-condition-differential-analysis.md), and the
[public API index](functions.md).
