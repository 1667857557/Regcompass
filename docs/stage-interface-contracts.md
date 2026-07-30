# Stage input-output contracts

RegCompass connects stages only when classes, workflow settings, GEM provenance,
analysis mode, metacell construction provenance and unit order agree. Equations
are in [Tutorial 3](tutorial-03-mathematical-model.md).

## Stage 1: Pando evidence

Class: `regcompass_grn_step`.

```r
step1$grn_result$analysis_mode
step1$grn_result$condition_coefficients_calculated
step1$grn_result$condition_fit_status
step1$grn_result$tf_peak_gene_condition
```

### `analysis_mode = "condition_grn"`

Selected when the effective condition column has at least two levels. The
canonical `pando_condition_grn_fit` contract contains:

```r
fit$coefficient_estimable_mask
fit$projectable_structural_zero_mask
fit$projection_support_mask
fit$projection_condition_full_oof
fit$projection_common_oof
fit$projection_global_common_oof
```

`coefficient_estimable_mask` records whether a coefficient can be fitted.
`projectable_structural_zero_mask` records a shared candidate edge whose
contribution is fixed at zero in that condition. The two masks are mutually
exclusive, and their union is `projection_support_mask`.

`projection_condition_full_oof` is the primary projection. Jointly estimable
edges form `projection_common_oof`; their difference is the condition-unique
component.

### `analysis_mode = "standard_pando"`

Selected when the condition column is omitted, absent or single-level.
RegCompass calls `Pando::infer_grn()` independently within each broad cell type.
No `ConditionGRNFit`, condition coefficient, condition deviation or condition
contrast is calculated.

## Stage 2: independent cell-type graphs, joint conditions

Class: `regcompass_metacell_step`.

```r
step2$pooled$metacell_meta
step2$pooled$membership
step2$pooled$input_design
step2$metacell_object
```

RegCompass scales RNA and ATAC embedding blocks within each broad cell type using
all conditions and calls:

```r
SuperCell::SCimplify_by_graph_group_from_embedding(
  X = joint_embedding,
  cell.graph.group = cell_type,
  cell.split.condition = condition,
  gamma = gamma
)
```

The contract is:

```text
native_supercell_api = SCimplify_by_graph_group_from_embedding
graph_group_argument = cell.graph.group
condition_argument = cell.split.condition
graph_scope = one_independent_graph_per_cell_type
condition_scope = all_conditions_joint_within_cell_type_graph
membership_split_timing = after_joint_graph_clustering
embedding_scaling = within_celltype_joint_condition_equal_modality_blocks
temporary_combined_stratum = FALSE
```

Final metacells are pure for cell type and condition. No sample-derived grouping
or concatenated condition-by-cell-type field is created. RNA and ATAC counts are
aggregated from exact `membership(cell_id, metacell_id)`.

## Stage 3: biological meta-modules

Class: `regcompass_meta_module_step`.

Active Pando target genes form one supported set per effective condition and cell
type. Positive and negative coefficients both count as evidence. A reaction is
core only when one complete GPR branch is represented. Stage 3 creates the
reaction catalogue and does not run FASTCORE.

## Stage 4: regulatory Layer 1

Class: `regcompass_layer1_step`.

Both modes use cell-first TF RNA × peak ATAC projection followed by exact
SuperCell aggregation. Interactions are never reconstructed from metacell means.

Condition mode follows:

```text
condition-full outer-heldout projection (primary)
+ common-support outer-heldout component
+ condition-unique projection difference
→ metacell mean
→ cell-type latent RNA support
→ reliability × tanh(primary projection / shared scale)
→ bounded RNA-support odds modifier
→ GPR reaction expression
```

Required schema fields include:

```r
step4$reaction_expression_condition_full_oof
step4$reaction_expression_common_oof
step4$reaction_expression_rna_only
step4$gene_projection_condition_full_oof
step4$gene_projection_common_oof
step4$gene_projection_condition_unique_oof
step4$projection_provenance
```

`reaction_expression` is identical to
`reaction_expression_condition_full_oof`. A non-estimable edge side contributes
zero. A non-finite target modifier uses neutral `R = 0`, exactly recovering
RNA-only support.

The Stage 4 schema does not contain depth-matching, common-depth,
alpha-sensitivity, zero-support-sensitivity or link-saturation-propagation
fields.

## Stage 5: shared model and directional penalties

Class: `regcompass_layer2_step`.

For each medium, one shared union GEM and one global FASTCORE completion are
reused for every condition, metacell and evidence route.

```r
step5$penalty_condition_full_oof
step5$penalty_common_oof
step5$penalty_condition_unique_increment
step5$penalty_rna_only
step5$vmax
step5$model_cache_summary
step5$structural_model_contract
```

`penalty` is identical to `penalty_condition_full_oof`. The condition-unique
increment is the primary penalty minus the common-support penalty. All routes
share reaction order, bounds, target direction and `vmax`.

## Stage 6: final result

```r
result$analysis_mode
result$reaction_ranking
result$condition_summary
result$condition_contrast
result$common_support_component_summary
result$condition_unique_penalty_increment_summary
result$rna_only_control_summary
```

Primary rankings and condition statistics use condition-full OOF. For one
effective condition, `condition_contrast` is empty and no artificial second
condition is generated.

Public API: [functions.md](functions.md).
