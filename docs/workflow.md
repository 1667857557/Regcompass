# RegCompassR workflow

This page describes the six canonical stages and the optional targeted reaction
remapping pass. Equations are in
[Tutorial 3](tutorial-03-mathematical-model.md).

## Data flow

```text
paired RNA+ATAC cells
→ condition-aware or standard Pando
→ independent cell-type graphs with joint conditions
→ condition-pure metacells
→ supported metabolic genes and reaction catalogue
→ condition-full regulatory reaction support
→ shared medium-specific metabolic model
→ directional penalties and condition comparisons
→ optional direct database-linked targeted remapping
```

## Stage 1: `rc_regcompass_step_grn()`

Targets are GEM GPR genes present in the RNA assay. With at least two conditions,
Pando fits one shared candidate supergraph per broad cell type using
`candidate_screen = "motif_domain"`, equal condition weights and nested
outer-heldout projection. Otherwise original `Pando::infer_grn()` is used and No
condition coefficients are calculated.

Main condition-mode outputs:

```r
step1$grn_result$condition_grn_fits
step1$grn_result$condition_fit_status
step1$grn_result$tf_peak_gene_condition
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
reaction requires one complete GPR branch. The catalogue adds direct subsystem,
KEGG/Reactome and master-Rhea relations, then merges condition- and cell-type-
specific sets. Stage 3 does not run FASTCORE.

## Stage 4: `rc_regcompass_step_layer1()`

The primary regulatory input is condition-full OOF projection. Jointly estimable
edges form the common-support component; each non-estimable edge side contributes
zero. Single-cell scores are averaged by exact metacell membership before RNA
support modification and GPR aggregation.

```r
step4$gene_projection_condition_full_oof
step4$gene_projection_common_oof
step4$gene_projection_condition_unique_oof
step4$reaction_expression_condition_full_oof
step4$reaction_expression_common_oof
```

The Stage 4 schema contains no depth-matching, common-depth, alpha-sensitivity,
zero-support-sensitivity or link-saturation-propagation branches.

## Stage 5: `rc_regcompass_step_layer2()`

One global FASTCORE completion is run per medium, and the completed model is
reused for every condition, metacell and evidence route. Forward and reverse
directions are scored separately.

```r
step5$penalty_condition_full_oof
step5$penalty_common_oof
step5$penalty_condition_unique_increment
step5$penalty_rna_only
step5$vmax
step5$model_cache_summary
```

The condition-full penalty is primary. Lower normalized penalty indicates
stronger network-constrained support; it is not measured flux.

## Stage 6: `rc_regcompass_step_results()`

```r
result$reaction_ranking
result$condition_summary
result$condition_contrast
result$common_support_component_summary
result$condition_unique_penalty_increment_summary
```

Condition comparisons must fix reaction, direction, medium, broad cell type,
model, bounds and target-flux fraction. Metacell P values describe within-dataset
separation and are not donor-level inference.

## Optional targeted remapping: `rc_regcompass_step_target_union()`

After Stage 5, selected reaction anchors can be used to identify directly linked
non-core reactions sharing KEGG, Reactome or master-Rhea identifiers. The second
pass reuses the exact cached medium-specific union GEMs and does not rerun
FASTCORE or reconstruct a model.

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
an optional target-extension analysis, not one of the removed sensitivity or
comparability guardrails.

## Restart boundaries

| Earliest stage | Changes |
|---|---|
| Stage 1 | RNA/ATAC data, labels, genome, regions, motifs or Pando fitting |
| Stage 2 | reductions, dimensions, gamma, seed, thresholds or cells |
| Stage 3 | GPR rules or reaction annotations |
| Stage 4 | projection, RNA support or GPR aggregation |
| Stage 5 | medium, bounds, direction, omega, solver or model completion |
| Stage 6 | annotations or reporting filters |
| Targeted remapping | selected anchors or direct cross-reference target set only |

Tutorial: [targeted reaction remapping](tutorial-04-targeted-reaction-remapping.md).
Public API: [functions.md](functions.md).
