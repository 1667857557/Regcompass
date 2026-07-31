# RegCompassR public API

Main workflow tutorials: [one-shot](tutorial-01-quick-start.md),
[stepwise](tutorial-02-stepwise-audit.md),
[mathematical model](tutorial-03-mathematical-model.md),
[targeted reaction remapping](tutorial-04-targeted-reaction-remapping.md), and
[condition comparison](tutorial-05-condition-differential-analysis.md).

## Complete workflows

| Function | Purpose |
|---|---|
| `rc_run_regcompass()` | Run all six stages with automatic standard/condition-aware Pando routing. |
| `rc_run_regcompass_one_shot()` | Convenience wrapper around the complete workflow with species-aware plasma defaults. |

`condition_col` may be absent, single-level or multi-level. The selected route is
returned in `result$analysis_mode`: `standard_pando` or `condition_grn`.

## Restartable stages

| Stage | Function | Main output |
|---:|---|---|
| 1 | `rc_regcompass_step_grn()` | standard Pando networks or canonical `pando_condition_grn_fit` contracts |
| 2 | `rc_regcompass_step_metacells()` | cell-type-scoped joint-condition WNN graphs and condition-pure metacells |
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
uses `Pando::infer_grn()` and produces no condition coefficients.

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
  k.knn = 30L,
  seed = 12345L,
  min_cells_per_stratum = 20L,
  min_metacell_size = 1L,
  min_metacells_per_stratum = 1L
)
```

RegCompass calls `SuperCell::SCimplify_by_graph_group()` with broad cell type as
`cell.graph.group` and condition as `cell.split.condition`. Each broad cell type
gets one independent native RNA+ATAC WNN graph. All conditions within that cell
type jointly determine adaptive modality weights, neighbours and Walktrap
parent clusters; condition splits membership only after clustering. No sample
column or concatenated condition-by-cell-type stratum is used.

The canonical default is `gamma = 30`. Small condition-split metacells are
retained and marked by `low_power_metacell` instead of being silently merged or
removed.

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

## Optional targeted reaction remapping

| Function | Purpose |
|---|---|
| `rc_regcompass_step_target_union()` | Score direct KEGG/Reactome/master-Rhea-linked non-core reactions in the exact cached Stage 5 union GEMs. |

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

The function does not rebuild the model or rerun FASTCORE. In condition mode it
uses `step4$reaction_expression`, the canonical alias of
`reaction_expression_condition_full_oof`, so targeted reactions remain on the
same primary evidence scale as the original Stage 5 targets.

## GEM and medium

| Function | Purpose |
|---|---|
| `rc_prepare_gem()` | Load and validate a pinned Human-GEM or Mouse-GEM. |
| `rc_make_medium_scenarios()` | Build authoritative-journal plasma or culture-challenge scenarios and user-defined exchange bounds. |
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

Human composition uses HPLM from *Cell* 2017 and updated HPLM from *Cell
Metabolism* 2021. Plasmax from *Science Advances* 2019 is an independent
validation source and is not numerically averaged with HPLM. Mouse composition
uses the package's documented mouse-plasma evidence table.

All five challenges use `medium_background_id = authoritative_HPLM_2017_2021`
and override only the named nutrient. Their construction is
`authoritative_HPLM_background_plus_named_nutrient_override`. Output provenance
includes `background_reference_doi`, `background_validation_reference_doi`,
`challenge_reference_doi`, and `scenario_construction`.

Users may supply reaction-level `custom_medium` or metabolite-level
`custom_metabolites` with `scenario = "custom"` or `scenario = NULL`. Technical
GEM boundary modes are not biological scenarios. The medium policy is
`authoritative_journal_composition_with_explicit_overrides`.

## Condition analysis

| Function | Purpose |
|---|---|
| `rc_test_condition_reactions()` | Compare fixed reaction-direction targets. |
| `rc_report_condition_directions()` | Summarize forward and reverse targets. |
| `rc_plot_condition_reaction()` | Plot one reaction direction across conditions. |
| `rc_select_gene_reactions()` | Select scored reactions by metabolic gene. |

For one condition, `result$condition_contrast` is empty. Metacell tests are
within-dataset inference, not biological-replicate inference.
