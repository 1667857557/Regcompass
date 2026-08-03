# RegCompassR workflow

This page describes the six canonical stages and the optional targeted reaction
remapping pass. Equations are in
[Tutorial 3](tutorial-03-mathematical-model.md).

## Data flow

```text
paired RNA+ATAC cells
→ condition-aware fixed-dictionary or standard Pando
→ independent cell-type graphs with joint conditions
→ condition-pure metacells
→ condition-specific biological meta-modules
→ union of conditions within each cell type only
→ cell-type-specific regulatory reaction support
→ one union GEM and independent FASTCORE per cell type × medium
→ directional penalties and within-cell-type condition comparisons
→ optional cell-type-scoped targeted remapping
```

## Stage 1: `rc_regcompass_step_grn()`

Targets are GEM GPR genes present in the RNA assay. With at least two retained
conditions, Pando performs pooled and per-condition biological candidate
discovery within each broad cell type. Exact `(TF, peak, target)` triples are
unioned into one frozen dictionary, and every condition is fit with the same
unscaled Gaussian identity interaction model. Otherwise original
`Pando::infer_grn()` is used and no condition coefficients are calculated.

For multi-condition fits, RegCompass accepts only estimable coefficients with
within-condition BH `padj < 0.05`. The pooled stage supports candidate recall and
does not rescale condition coefficients.

Main condition-mode outputs:

```r
step1$grn_result$condition_grn_fits
step1$grn_result$condition_fit_status
step1$grn_result$tf_peak_gene_condition_effect_all
step1$grn_result$tf_peak_gene_condition_effect
```

## Stage 2: `rc_regcompass_step_metacells()`

RNA and ATAC embeddings are standardized within each cell type across all
conditions. SuperCell builds one independent graph per cell type and applies
condition after graph clustering to create condition-pure metacells. Exact
membership is used to aggregate RNA and ATAC counts.

```r
step2$pooled$membership
step2$pooled$metacell_meta
step2$pooled$input_design
```

## Stage 3: `rc_regcompass_step_meta_modules()`

Active condition coefficients identify supported metabolic genes. A core
reaction requires one complete GPR branch. For every condition-by-cell-type
group, the biological catalogue adds direct subsystem, KEGG/Reactome and
master-Rhea relations.

Condition-specific catalogues are then unioned **within the same cell type
only**. Different cell types remain in separate `cell_type_catalogues`; Stage 3
never creates a GEM and never runs FASTCORE.

```r
step3$merged_modules$cell_type_catalogues
step3$merged_modules$merged_core_reactions
step3$merged_modules$merged_reaction_membership
```

The merged tables retain the workflow cell-type column and record
`merge_scope = "cell_type"` and `cross_celltype_merge = FALSE`.

## Stage 4: `rc_regcompass_step_layer1()`

For each condition, significant fixed-dictionary edge effects are projected to
paired cells as `penalty_effect × TF_RNA × peak_ATAC`, summed by target gene and
averaged using exact metacell membership. RNA support and GPR rules then produce
reaction-level evidence.

```r
step4$gene_projection_condition_full_oof
step4$gene_projection_common_oof
step4$gene_projection_condition_unique_oof
step4$reaction_expression_condition_full_oof
step4$reaction_expression_common_oof
```

These historical field names are compatibility aliases. The current estimator
is not OOF: `condition_full_oof` is the primary fixed-dictionary route,
`common_oof` aliases the primary route, and `condition_unique_oof` is a zero
compatibility matrix. Targets without significant estimable regulatory edges
use the neutral RNA-only fallback.

## Stage 5: `rc_regcompass_step_layer2()`

For `model_mode = "meta_module_gem"`, Stage 5 builds one union GEM for every
`cell_type × medium_scenario` pair. Conditions of that cell type contribute to
its biological reaction union; reactions from other cell types do not.
FASTCORE runs independently within each cell-type/medium model.

Each model is reused only for metacells whose cell type matches the model.
Forward and reverse directions are scored separately. Directional `vmax` is
computed once per cell-type model and target direction, then reused across
conditions and metacells of that cell type.

```r
step5$penalty_condition_full_oof
step5$penalty_common_oof
step5$penalty_condition_unique_increment
step5$penalty_rna_only
step5$vmax
step5$model_cache_summary
step5$structural_model_contract
```

The cache summary carries `cell_type`, medium, model file, checksum, reaction
counts and cell-type FASTCORE support counts. The condition-specific penalty is
primary. Lower normalized penalty indicates stronger network-constrained
support; it is not measured flux.

The optional `full_gem` mode uses the full reference GEM and remains separate
from cell-type union-GEM construction.

## Stage 6: `rc_regcompass_step_results()`

```r
result$reaction_ranking
result$condition_summary
result$condition_contrast
result$common_support_component_summary
result$condition_unique_penalty_increment_summary
```

Condition comparisons fix reaction, direction, medium and broad cell type.
Compared conditions share the same cell-type/medium model, bounds and target
`vmax`. Cross-cell-type rows are excluded rather than interpreted as missing
observations on a global model. Metacell P values describe within-dataset
separation and are not donor-level inference.

## Optional targeted remapping: `rc_regcompass_step_target_union()`

After Stage 5, selected reaction anchors can identify directly linked non-core
reactions sharing KEGG, Reactome or master-Rhea identifiers. The second pass
reuses the exact cached union GEMs for the corresponding cell type and medium.
Availability is intersected across media within each cell type; models from
different cell types are never merged. FASTCORE is not rerun.

```r
targeted <- rc_regcompass_step_target_union(
  layer1 = step4,
  meta_modules = step3,
  layer2 = step5,
  gem = gem,
  outdir = "RegCompass_targeted",
  core_reaction_ids = c("MAR04381", "MAR04379"),
  layer2_args = list(target_direction = "both", solver = "highs")
)
```

In condition mode, `step4$reaction_expression` is the canonical alias of
`reaction_expression_condition_full_oof`; targeted reactions therefore use the
same primary regulatory evidence route as the original Stage 5 scoring. This is
an optional target-extension analysis, not a sensitivity or comparability
branch.

## Restart boundaries

| Earliest stage | Changes |
|---|---|
| Stage 1 | RNA/ATAC data, labels, genome, regions, motifs or Pando fitting |
| Stage 2 | reductions, dimensions, gamma, seed, thresholds or cells |
| Stage 3 | GPR rules, reaction annotations or cell-type module membership |
| Stage 4 | projection, RNA support or GPR aggregation |
| Stage 5 | medium, bounds, direction, omega, solver or cell-type FASTCORE completion |
| Stage 6 | annotations or reporting filters |
| Targeted remapping | selected anchors or direct cross-reference target set only |

Tutorial: [targeted reaction remapping](tutorial-04-targeted-reaction-remapping.md).
Public API: [functions.md](functions.md).
