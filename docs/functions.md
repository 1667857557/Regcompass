# Public functions in RegCompassR 1.8.9

## Setup and complete runs

- `rc_prepare_gem()`, `rc_prepare_human2_gem()`, `rc_prepare_mouse_gem()`: load bundled pinned GEMs or explicitly rebuild an upstream release.
- `rc_bundled_gem_manifest()`: inspect bundled model provenance.
- `rc_make_medium_scenarios()`: create the shared medium table.
- `rc_run_regcompass()` and `rc_run_regcompass_one_shot()`: execute the six-stage shared-background condition-sub-GRN workflow.
- `rc_parallel_config()`: inspect operating-system backend resolution.

Canonical GRN mode:

```r
grn_mode = "multitask_shared_backbone"
```

Legacy reproducibility mode:

```r
grn_mode = "legacy_condition_pando"
```

The canonical complete runner and Stage 1/2 functions require `condition_col` and `celltype_col` only. They do not expose a biological-sample column.

## Inspectable stages

- `rc_regcompass_step_grn()`: build one validated Pando structural TF–peak–target universe per cell type, apply a condition-aware observability filter, estimate global and condition coefficients with condition-balanced elastic net, and calculate bootstrap reproducibility.
- `rc_regcompass_step_metacells()`: build condition-stratified, cell-type-labelled SuperCell2 metacells.
- `rc_regcompass_step_meta_modules()`: map active condition sub-GRN targets to complete-GPR core reactions and one ordered subsystem/KEGG–Reactome/master-Rhea expansion.
- `rc_regcompass_step_layer1()`: combine RNA support with an ATAC-only projection of the fitted condition sub-GRN and aggregate through GPR rules.
- `rc_regcompass_step_layer2()`: build one medium-specific union GEM, apply one global FASTCORE completion, and run directional COMPASS-like LP scoring.
- `rc_regcompass_step_results()`: assemble compact rankings, contrasts, regulatory evidence, complete-GPR cores, and provenance.
- `rc_regcompass_step_target_union()`: map additional targets and score them in the exact cached Stage 5 model.

## Stage 1 structural design controls

```r
pando_args = list(
  min_cells = 100L,
  pando_design_args = list(
    peak_to_gene_method = "Signac",
    upstream = 100000,
    downstream = 0,
    extend = 1000000,
    only_tss = FALSE,
    min_tf_detection = 0,
    min_peak_detection = 0,
    min_target_detection = 0,
    max_edges_per_target = Inf
  )
)
```

`Pando::prepare_grn_design()` creates a version-2 condition-agnostic dictionary. Exact predictors are deduplicated by TF, measured ATAC feature, and target. `supporting_regions` retains multiple regulatory regions mapping to one peak. RegCompass validates the design fingerprint and feature mapping before fitting.

The pooled Pando detection thresholds remain zero so a condition-restricted TF, peak, or target is not removed before the shared model is formed. A finite `max_edges_per_target` is rejected because Pando candidate order is deterministic but not an evidence ranking.

## Stage 1 multitask, observability, and bootstrap controls

```r
multitask_args = list(
  alpha = 0.5,
  global_penalty_factor = 1,
  deviation_penalty_factor = 1,
  lambda_rule = "lambda.1se",
  nfolds = 5L,
  n_bootstrap = 100L,
  min_selection_frequency = 0.7,
  min_sign_stability = 0.8,
  min_abs_effect = 0,
  min_cv_rsq = 0,
  min_bootstrap_success_fraction = 0.8,
  candidate_screen_threshold = 0,
  max_edges_per_target = Inf,
  min_detected_cells_per_condition = 10L,
  min_detection_fraction_per_condition = 0.01,
  seed = 12345L
)
```

For condition \(c\), RegCompass uses

\[
m_c=\min\left(n_c,\max\left(10,\left\lceil0.01n_c\right\rceil\right)\right).
\]

An edge enters the shared model universe only when the non-zero TF-RNA × peak-ATAC predictor and target RNA each occur in at least \(m_c\) cells of one or more conditions. TF and peak must be non-zero in the same cells. This filter uses detection only and does not use target correlation or fitted effect size.

`alpha = 0.5` retains a lasso component for sparse selection and a ridge component for correlated TF–peak predictors. Global and zero-sum condition-deviation coordinates use equal explicit penalty factors by default. `deviation_penalty_factor > 1` is an optional sensitivity prior favoring a more conserved shared backbone.

Five-fold cross-validation estimates condition centres and edge scales from each training fold. An active target requires strictly positive out-of-fold `cv_rsq`; a positive `min_cv_rsq` adds a stronger floor.

`n_bootstrap = 100` controls full-size nonparametric bootstrap fits. Every condition is sampled with replacement at its original cell count. Bootstrap targets and predictors are re-centred within condition before fitting at the full-data lambda and scale. Edge tables report selection-frequency Monte Carlo standard errors and Wilson 95% intervals.

`min_sign_stability = 0.8` corresponds to at least 90% agreement on one sign among selected bootstrap fits because \(\rho=|2q-1|\).

Principal Stage 1 outputs:

```r
step1$grn_result$tf_peak_gene_candidates
step1$grn_result$tf_peak_gene_global
step1$grn_result$tf_peak_gene_condition_all
step1$grn_result$tf_peak_gene_significant
step1$grn_result$condition_target_genes
step1$grn_result$target_model_diagnostics
step1$grn_result$stability_diagnostics
step1$grn_result$group_status
```

The candidate table distinguishes:

```text
edge_universe_id       = complete Pando structural universe
model_edge_universe_id = observability-filtered shared model universe
model_observable       = whether a structural edge entered fitting
```

## Stage 2 geometry

```r
metacell_args = list(
  rna_reduction = "pca",
  rna_dims = 1:30,
  atac_reduction = "lsi",
  atac_dims = 2:30,
  gamma = 30L,
  seed = 12345L,
  min_cells_per_stratum = 100L,
  min_metacell_size = 20L,
  min_metacells_per_stratum = 2L
)
```

Condition is the RegCompass stratum. Cell type is passed to SuperCell2 as the exact label. An existing Harmony reduction may replace PCA. Changing cells, assays, reductions, dimensions, seed, gamma, or thresholds invalidates the Stage 2 cache.

## Stage 3 condition and merged catalogues

```r
step3$condition_modules$supported_metabolic_genes
step3$condition_modules$core_gene_reaction
step3$condition_modules$reaction_membership
step3$merged_modules$merged_core_reactions
step3$merged_modules$merged_reaction_membership
```

`group_id` is the `condition × cell type` analysis identifier. A reaction becomes core only when one complete GPR branch is contained in the condition target set. The merged object is a biological reaction catalogue, not a GEM.

## Stage 4 integrated evidence

```r
layer1_args = list(
  regulatory_alpha = 1,
  gpr_and_method = "min",
  gene_half_saturation = 1
)
```

`gpr_and_method` accepts `min`, `median`, or `mean`. Isozyme OR branches are additive. A gene without an active condition edge has a zero regulatory modifier and exactly returns its RNA-only support.

## Stage 5 union-GEM controls

```r
layer2_args = list(
  target_direction = "both",
  solver = "highs",
  model_params = list(
    completion_time_limit = 600,
    fastcore_epsilon = 1e-4,
    max_support_reactions = 2000,
    strict = TRUE
  )
)
```

`completion_time_limit` applies only to union-GEM construction. Directional scoring LPs have no time-limit argument. Every condition and metacell under one medium uses the same cached model file, reaction IDs, stoichiometry, and bounds.

## Interpretation and plotting

- `rc_build_reaction_annotations()` and `rc_attach_reaction_annotations()`: attach reaction names, formulas, GPRs, and database cross-references.
- `rc_test_condition_reactions()`: descriptive same-reaction, same-direction, same-medium comparisons within cell type.
- `rc_report_condition_directions()`: retain forward/reverse targets and report direction-aware support without claiming net flux.
- `rc_select_gene_reactions()`: select scored reactions by GPR gene.
- `rc_plot_condition_reaction()` and `rc_plot_condition_gene_reactions()`: condition plots with reaction and evidence annotations.

Bootstrap measures selection reproducibility under cell resampling; it does not create biological-replicate inference. GRN fitting and RNA support use the same paired multiome dataset. Condition centring and ATAC-only projection reduce direct duplicate weighting but do not create independent validation evidence. Metacell-level tests remain descriptive pseudo-observation analyses.

## Documentation

- [Pando and multitask GRN parameter policy](grn-parameter-policy.md)
- [Mathematics and object contracts](multitask-shared-grn.md)
- [Quick start](tutorial-01-quick-start.md)
- [Stepwise audit](tutorial-02-stepwise-audit.md)
- [Restart and diagnostics](tutorial-03-advanced-restart.md)
- [Targeted reaction remapping](tutorial-04-targeted-reaction-remapping.md)
- [Condition differential analysis](tutorial-05-condition-differential-analysis.md)
- [Stage contracts](stage-interface-contracts.md)
- [Portable execution](portable-execution.md)
