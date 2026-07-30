# RegCompassR public API

Main workflow tutorials: [one-shot](tutorial-01-quick-start.md),
[stepwise](tutorial-02-stepwise-audit.md),
[mathematical model](tutorial-03-mathematical-model.md), and
[condition comparison](tutorial-04-condition-differential-analysis.md).

## Complete workflows

| Function | Purpose |
|---|---|
| `rc_run_regcompass()` | Run all six stages with automatic standard/condition-aware Pando routing. |
| `rc_run_regcompass_one_shot()` | Convenience wrapper around the complete workflow. |

`condition_col` may be absent, single-level or multi-level. The selected route is
returned in `result$analysis_mode`: `standard_pando` or `condition_grn`.

## Restartable stages

| Stage | Function | Main output |
|---:|---|---|
| 1 | `rc_regcompass_step_grn()` | standard Pando networks or canonical `pando_condition_grn_fit` contracts |
| 2 | `rc_regcompass_step_metacells()` | independent cell-type graphs and condition-pure metacells |
| 3 | `rc_regcompass_step_meta_modules()` | supported genes, complete-GPR cores and reaction catalogue |
| 4 | `rc_regcompass_step_layer1()` | condition-full, common-support and RNA reaction expression |
| 5 | `rc_regcompass_step_layer2()` | shared structural model and directional penalties |
| 6 | `rc_regcompass_step_results()` | annotations, primary rankings and decomposition outputs |

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

In condition mode, `candidate_screen = "motif_domain"`,
`condition_weight = "equal"`, and `scale = TRUE` are required. Standard mode
uses `Pando::infer_grn()` and produces No condition coefficients.

Condition-mode projection contracts expose:

```r
fit$coefficient_estimable_mask
fit$projectable_structural_zero_mask
fit$projection_support_mask
fit$projection_condition_full_oof
fit$projection_common_oof
```

## Stage 2 graph contract

```r
metacell_args = list(
  rna_reduction = "pca",
  rna_dims = 1:30,
  atac_reduction = "lsi",
  atac_dims = 2:30,
  gamma = 30L,
  seed = 12345L
)
```

RegCompass calls `SuperCell::SCimplify_by_graph_group_from_embedding()` with
`cell.graph.group = broad cell type` and
`cell.split.condition = condition`. This produces
`one_independent_graph_per_cell_type` while preserving
`all_conditions_joint_within_cell_type_graph`; `temporary_combined_stratum = FALSE`.

## Stage 4 condition-full support

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

- `gene_projection_condition_full_oof` is primary;
- `gene_projection_common_oof` is the jointly estimable component;
- `gene_projection_condition_unique_oof` is their difference;
- a non-estimable edge side contributes a projectable structural zero;
- a non-finite target modifier uses neutral `R = 0` and equals RNA-only support;
- GPR AND uses `min` by default and OR isozyme branches are additive.

## Stage 5 penalty outputs

```r
step5$penalty_condition_full_oof
step5$penalty_common_oof
step5$penalty_condition_unique_increment
step5$penalty_rna_only
```

All four matrices share the same medium-specific GEM, bounds, reaction order,
target direction and `vmax`. The first matrix is the primary penalty.

The schema does not contain depth-matching, common-depth, alpha-sensitivity,
zero-support-sensitivity or link-saturation-propagation outputs.

## GEM and medium

| Function | Purpose |
|---|---|
| `rc_prepare_gem()` | Load and validate a pinned Human-GEM or Mouse-GEM. |
| `rc_make_medium_scenarios()` | Build preset or custom exchange bounds. |
| `rc_build_reaction_annotations()` | Build reaction names, formulas, GPRs and cross-references. |
| `rc_attach_reaction_annotations()` | Attach annotations to an existing result. |

## Condition analysis

| Function | Purpose |
|---|---|
| `rc_test_condition_reactions()` | Compare fixed reaction-direction targets. |
| `rc_report_condition_directions()` | Summarize forward and reverse targets. |
| `rc_plot_condition_reaction()` | Plot one reaction direction across conditions. |
| `rc_select_gene_reactions()` | Select scored reactions by metabolic gene. |

For one condition, `result$condition_contrast` is empty. Metacell tests are
within-dataset inference, not biological-replicate inference.
