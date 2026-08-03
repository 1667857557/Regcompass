# RegCompassR public API

Main tutorials: [one-shot](tutorial-01-quick-start.md),
[stepwise](tutorial-02-stepwise-audit.md),
[condition-GRN contract](condition-comparable-grn.md),
[targeted reaction remapping](tutorial-04-targeted-reaction-remapping.md), and
[condition comparison](tutorial-05-condition-differential-analysis.md).

## Complete workflows

| Function | Purpose |
|---|---|
| `rc_run_regcompass()` | Run all six stages with automatic standard/common-dictionary Pando routing. |
| `rc_run_regcompass_one_shot()` | Convenience wrapper with species-aware plasma defaults. |

`condition_col` may be absent, single-level or multi-level. The selected route is
returned in `result$analysis_mode`: `standard_pando` or `condition_grn`.

## Restartable stages

| Stage | Function | Main output |
|---:|---|---|
| 1 | `rc_regcompass_step_grn()` | standard Pando networks or `pando_condition_grn_common_dictionary_v1` contracts |
| 2 | `rc_regcompass_step_metacells()` | cell-type-scoped joint-condition WNN graphs and condition-pure metacells |
| 3 | `rc_regcompass_step_meta_modules()` | supported genes, complete-GPR cores and reaction catalogue |
| 4 | `rc_regcompass_step_layer1()` | paired-cell regulatory projection and reaction expression |
| 5 | `rc_regcompass_step_layer2()` | shared structural model and directional penalties |
| 6 | `rc_regcompass_step_results()` | annotations, rankings and condition contrasts |

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
    min_cells = 300L,
    min_abs_estimate = 0,
    min_model_rsq = 0.1,
    pando_infer_args = list(
      tf_cor = 0.1,
      peak_cor = 0,
      adjust_method = "BH",
      padj_threshold = 0.05,
      rank_action = "mark",
      min_residual_df = 1L
    )
  )
)
```

In condition mode, Pando performs global plus each-condition candidate discovery,
unions exact `(target, TF, region)` triples, freezes the dictionary, and fits the
same unscaled Gaussian identity interaction model in every condition. Condition
effects are the fitted coefficients; no pooled-coefficient calibration is used.

The following condition controls are retired and rejected:

```text
candidate_screen
condition_mix
condition_weight
alpha
nlambda / lambda / lambda_min_ratio
outer_nfolds / inner_nfolds
lambda_selection
scale
engine_control
comparison_conditions
```

Standard mode uses original `Pando::infer_grn()` independently within each broad
cell type and creates no condition coefficients.

Condition fit contracts expose:

```r
fit$edge_dictionary
fit$coefficients
fit$fit
fit$condition_cell_ids
fit$padj_threshold
fit$projection_effect_column
fit$projection_policy
```

The coefficient table retains `estimate`, `std_err`, `statistic`, `pval`, `padj`,
`significant`, `penalty_effect`, `estimable`, `zero_variance` and `aliased`.
`penalty_effect` equals `estimate` only for `padj < 0.05`.

## Stage 2 graph contract

```r
metacell_args = list(
  rna_reduction = "pca",
  rna_dims = 1:30,
  atac_reduction = "lsi",
  atac_dims = 2:30,
  gamma = 30L,
  k.knn = 30L,
  seed = 12345L
)
```

RegCompass calls `SuperCell::SCimplify_by_graph_group()` with broad cell type as
`cell.graph.group` and condition as `cell.split.condition`. Each broad cell type
gets one independent RNA+ATAC WNN graph. Conditions jointly determine that graph
and split membership only after clustering.

## Stage 4 regulatory support

```r
step4 <- rc_regcompass_step_layer1(
  grn,
  metacells,
  meta_modules,
  gem,
  outdir,
  projection_component = "condition",
  regulatory_alpha = 1,
  gpr_and_method = "min"
)
```

For multiple conditions, Pando reconstructs paired-cell `TF RNA × peak ATAC`,
applies `penalty_effect`, sums by target, and RegCompass aggregates by exact
SuperCell membership. No metacell-level coefficient fitting or TF×ATAC
reconstruction is performed.

Historical fields containing `_oof`, `common`, or `condition_unique` are retained
for API compatibility. The primary and common fields are aliases of the current
BH-filtered fixed-dictionary full-fit projection; the condition-unique
compatibility decomposition is zero.

A non-finite target modifier uses neutral `R = 0` and therefore RNA-only support.
GPR AND uses `min` by default and OR isozyme branches are additive.

## Stage 5 penalty outputs

```r
step5$penalty
step5$penalty_condition_full_oof
step5$penalty_common_oof
step5$penalty_condition_unique_increment
step5$penalty_rna_only
```

All matrices share the same medium-specific GEM, bounds, reaction order, target
direction and `vmax`. `penalty` is primary. Historical full/common fields are
compatibility aliases and the condition-unique increment is zero under the
current condition-GRN design.

## Optional targeted reaction remapping

| Function | Purpose |
|---|---|
| `rc_regcompass_step_target_union()` | Score direct KEGG/Reactome/master-Rhea-linked non-core reactions in the cached Stage 5 union GEMs. |

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

## GEM and medium

| Function | Purpose |
|---|---|
| `rc_prepare_gem()` | Load and validate a pinned Human-GEM or Mouse-GEM. |
| `rc_make_medium_scenarios()` | Build plasma, culture-challenge or custom exchange bounds. |
| `rc_build_reaction_annotations()` | Build reaction names, formulas, GPRs and cross-references. |
| `rc_attach_reaction_annotations()` | Attach annotations to an existing result. |

Supported biological medium identifiers:

```text
normal_human_plasma
mouse_plasma
high_glucose
low_glucose
high_lactate
low_lactate
low_glutamine
custom
```

## Condition analysis

| Function | Purpose |
|---|---|
| `rc_test_condition_reactions()` | Compare fixed reaction-direction targets. |
| `rc_report_condition_directions()` | Summarize forward and reverse targets. |
| `rc_plot_condition_reaction()` | Plot one reaction direction across conditions. |
| `rc_select_gene_reactions()` | Select scored reactions by metabolic gene. |

For one condition, `result$condition_contrast` is empty. Metacell tests are
within-dataset inference, not biological-replicate inference.
