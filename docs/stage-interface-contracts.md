# Stage input-output contracts

RegCompass connects stages only when classes, workflow settings, GEM provenance,
analysis mode, metacell construction provenance, cell-type scope and unit order
agree. Equations are in [Tutorial 3](tutorial-03-mathematical-model.md).

## Stage 1: Pando evidence

Class: `regcompass_grn_step`.

```r
step1$grn_result$analysis_mode
step1$grn_result$condition_coefficients_calculated
step1$grn_result$condition_fit_status
step1$grn_result$condition_grn_fits
step1$grn_result$tf_peak_gene_condition_effect_all
step1$grn_result$tf_peak_gene_condition_effect
```

### `analysis_mode = "condition_grn"`

Selected when the effective condition column has at least two retained levels.
For each broad cell type, the canonical Pando contract contains:

```r
fit$edge_dictionary
fit$coefficients
fit$condition_levels
fit$condition_cell_ids
fit$fit_engine
fit$coefficient_scale
fit$interaction
fit$projection_effect_column
```

The contract requires:

```text
schema = pando_condition_grn_common_dictionary_v1
fit_engine = two-stage exact-edge-union fixed-dictionary GLM
coefficient_scale = raw unscaled condition coefficient
interaction = TF:peak
```

Candidate discovery is performed in the pooled cell type and separately in each
condition. The exact observed `(TF, peak, target)` triples are unioned into one
frozen dictionary; no Cartesian edge expansion is permitted. Every condition
must contain every dictionary edge exactly once in the complete coefficient
table.

Within each condition, RegCompass recomputes BH values and requires:

```text
significant = estimable AND padj < 0.05
penalty_effect = estimate, when significant
penalty_effect = 0, otherwise
```

A non-estimable coefficient remains `NA` in the complete table. It is not a
structural zero and is excluded from regulatory evidence. Pooled coefficients
are not used to rescale condition coefficients.

### `analysis_mode = "standard_pando"`

Selected when the condition column is omitted, absent or single-level.
RegCompass calls `Pando::infer_grn()` independently within each broad cell type.
No condition coefficient, condition deviation or condition contrast is
calculated.

## Stage 2: cell-type-scoped WNN graphs, joint conditions

Class: `regcompass_metacell_step`.

```r
step2$pooled$metacell_meta
step2$pooled$membership
step2$pooled$input_design
step2$metacell_object
```

RegCompass calls the canonical SuperCell grouped builder once on the supplied
paired RNA+ATAC Seurat object:

```r
SuperCell::SCimplify_by_graph_group(
  seurat = object,
  cell.graph.group = cell_type,
  cell.split.condition = condition,
  assay = c(rna_assay, atac_assay),
  reduction = list(rna_reduction, atac_reduction),
  dims = list(rna_dims, atac_dims),
  gamma = 30,
  k.knn = 30
)
```

The contract is:

```text
native_supercell_api = SCimplify_by_graph_group
graph_group_argument = cell.graph.group
condition_argument = cell.split.condition
graph_method = multimodal_WNN
graph_scope = one_independent_WNN_graph_per_cell_type
condition_scope = all_conditions_joint_within_cell_type_graph
membership_split_timing = after_joint_WNN_graph_clustering
modality_weighting = adaptive_WNN_within_cell_type
temporary_combined_stratum = FALSE
```

For each broad cell type, all conditions jointly determine native RNA+ATAC WNN
modality weights, neighbours and Walktrap parent clusters. Condition is applied
only after clustering to obtain condition-pure final memberships. No sample
column or concatenated condition-by-cell-type stratum is used.

The canonical `gamma` default is 30. Final RNA and ATAC counts are aggregated
from the exact `membership(cell_id, metacell_id)` by a membership-mode
`SCimplify_for_Seurat()` call; the graph is not rebuilt during aggregation.
Small final metacells are retained and marked by `low_power_metacell` rather than
silently merged or removed.

## Stage 3: biological meta-modules

Class: `regcompass_meta_module_step`.

Active Pando target genes form one supported set per effective condition and cell
type. Positive and negative significant coefficients both count as evidence. A
reaction is core only when one complete GPR branch is represented. The catalogue
then adds the configured biological subsystem and direct database-equivalence
relations. Stage 3 does not run FASTCORE and does not construct a GEM.

Condition-specific catalogues are merged only within the same cell type. The
contract requires:

```r
step3$merged_modules$cell_type_catalogues
step3$merged_modules$merged_core_reactions
step3$merged_modules$merged_reaction_membership
```

```text
merge_scope = cell_type
cross_celltype_merge = FALSE
is_gem = FALSE
fastcore_applied = FALSE
```

Every merged core and membership row retains the workflow cell-type column.
Different cell types may contain different core and biological reaction sets.

## Stage 4: regulatory Layer 1

Class: `regcompass_layer1_step`.

Both condition-aware and standard routes use cell-first TF RNA × peak ATAC
projection followed by exact SuperCell aggregation. Interactions are never
reconstructed from metacell means.

Condition mode follows:

```text
significant estimable fixed-dictionary edges for one condition
→ penalty_effect × TF_RNA × peak_ATAC per paired cell
→ sum by target gene
→ exact metacell mean
→ cell-type latent RNA support
→ reliability × tanh(primary projection / shared within-cell-type scale)
→ bounded RNA-support odds modifier
→ GPR reaction expression
```

Required compatibility fields include:

```r
step4$reaction_expression_condition_full_oof
step4$reaction_expression_common_oof
step4$reaction_expression_rna_only
step4$gene_projection_condition_full_oof
step4$gene_projection_common_oof
step4$gene_projection_condition_unique_oof
step4$projection_provenance
```

The names containing `_oof`, `common` and `condition_unique` are retained for
schema compatibility. The current estimator is not OOF:

```text
condition_full_oof = primary fixed-dictionary projection
common_oof = compatibility alias of primary
condition_unique_oof = zero compatibility matrix
```

`reaction_expression` is identical to the primary fixed-dictionary route. A
target without a significant estimable regulatory edge uses neutral regulatory
modification, exactly recovering RNA-only support.

## Stage 5: cell-type structural models and directional penalties

Class: `regcompass_layer2_step`.

For `model_mode = "meta_module_gem"`, the structural key is:

```text
cell_type × medium_scenario
```

For every key, Stage 5:

1. takes only the merged biological reactions for that cell type;
2. applies the medium-specific bounds;
3. runs FASTCORE independently for that cell-type model;
4. writes a distinct model file and checksum;
5. computes directional `vmax` once per target direction;
6. reuses the model only for metacells whose cell type matches.

The contract requires:

```text
shared_across_conditions = TRUE
shared_across_cell_types = FALSE
structural_scope = cell_type_x_medium
completion_stage = celltype_specific_fastcore_after_condition_module_union
```

```r
step5$penalty_condition_full_oof
step5$penalty_common_oof
step5$penalty_condition_unique_increment
step5$penalty_rna_only
step5$vmax
step5$model_cache_summary
step5$structural_model_contract
```

`penalty` is identical to the primary fixed-dictionary penalty. The common field
is a compatibility alias and the condition-unique increment is a zero
compatibility decomposition. All evidence routes for a cell type share its exact
reaction order, bounds, target direction and `vmax`; no route shares a union GEM
with another cell type.

The optional `full_gem` mode uses the complete reference GEM and is dispatched to
a separate full-GEM engine. It does not construct a cross-cell-type union GEM.

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

Primary rankings and condition statistics use the fixed-dictionary condition
route. Comparisons fix cell type, reaction, direction and medium; rows from
another cell type are excluded. For one effective condition,
`condition_contrast` is empty and no artificial second condition is generated.

Public API: [functions.md](functions.md).
