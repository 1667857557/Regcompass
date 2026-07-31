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
projected to paired cells and aggregated to metacells. No condition coefficient
or condition contrast is calculated.

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

After Stage 5, `rc_regcompass_step_target_union()` remains available as an
optional targeted-remapping pass. It reuses the cached Stage 5 union GEM and the
canonical Layer 1 `reaction_expression` route; in condition mode that route is
`reaction_expression_condition_full_oof`.

## SuperCell graph and purity inputs

Stage 2 calls `SuperCell::SCimplify_by_graph_group()` with:

```text
seurat                  = paired RNA+ATAC Seurat object
cell.graph.group        = broad cell type
cell.split.condition    = condition
gamma                   = 30 by default
k.knn                   = 30 by default
```

Each broad cell type receives one independent native multimodal WNN graph. All
conditions within that cell type jointly determine adaptive RNA/ATAC modality
weights, neighbours and Walktrap parent clusters. `cell.split.condition` is
applied only after clustering, so final metacells are condition-pure without
separate condition graphs.

No combined metadata stratum or sample-derived grouping is created. Provenance
records:

```text
native_supercell_api = SCimplify_by_graph_group
graph_scope = one_independent_WNN_graph_per_cell_type
condition_scope = all_conditions_joint_within_cell_type_graph
membership_split_timing = after_joint_WNN_graph_clustering
modality_weighting = adaptive_WNN_within_cell_type
temporary_combined_stratum = FALSE
```

Final RNA and ATAC counts are aggregated using the exact final membership.
Small condition-split metacells are retained and marked rather than silently
merged or removed.

## Restart boundaries

- motifs, regions, targets or Pando fitting changed: rerun Stage 1 onward;
- reductions, dimensions, gamma or SuperCell settings changed: rerun Stage 2 onward;
- GPR or catalogue annotations changed: rerun Stage 3 onward;
- projection, RNA support or GPR aggregation changed: rerun Stage 4 onward;
- medium, union GEM, FASTCORE or LP controls changed: rerun Stage 5 onward;
- targeted anchors or direct cross-reference targets changed: rerun only targeted remapping.

## Result behavior

Multiple conditions produce primary condition-full summaries and pairwise
contrasts. Common-support and RNA-only routes remain decomposition/control
outputs. A single effective condition produces rankings and summaries with an
empty condition contrast.

See [Tutorial 1](tutorial-01-quick-start.md),
[Tutorial 2](tutorial-02-stepwise-audit.md),
[Tutorial 3](tutorial-03-mathematical-model.md),
[Tutorial 4](tutorial-04-targeted-reaction-remapping.md),
[Tutorial 5](tutorial-05-condition-differential-analysis.md), and the
[public API index](functions.md).
