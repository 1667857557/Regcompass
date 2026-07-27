# Public functions in RegCompassR 1.8.10

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

Recommended metadata input:

```r
condition_col = "condition"
celltype_col = "cell_type"
sample_col = "sample_id"
```

`condition_col` and `celltype_col` are required. `sample_col` is optional and is used only for Stage 1 bootstrap. A valid sample column activates condition-stratified sample/donor cluster resampling. When the argument is omitted or names an absent column, RegCompass prints the exact fallback reason and uses condition-stratified cell resampling. Stage 2 does not expose `sample_col` and remains condition-only.

## Inspectable stages

- `rc_regcompass_step_grn()`: build one validated GREAT-domain Pando structural TF–peak–target universe per cell type, apply a condition-aware observability filter, estimate direct condition-specific coefficients, and calculate sample-aware bootstrap reproducibility.
- `rc_regcompass_step_metacells()`: build condition-stratified, cell-type-labelled SuperCell2 metacells.
- `rc_regcompass_step_meta_modules()`: map active condition sub-GRN targets to complete-GPR core reactions and one ordered subsystem/KEGG–Reactome/master-Rhea expansion.
- `rc_regcompass_step_layer1()`: combine RNA support with an ATAC-only projection of the fitted condition sub-GRN and aggregate through GPR rules.
- `rc_regcompass_step_layer2()`: build one medium-specific union GEM, apply one global FASTCORE completion, and run directional COMPASS-like LP scoring.
- `rc_regcompass_step_results()`: assemble compact rankings, contrasts, regulatory evidence, complete-GPR cores, bootstrap provenance, and structural provenance.
- `rc_regcompass_step_target_union()`: map additional targets and score them in the exact cached Stage 5 model.

## Stage 1 structural design controls

```r
pando_args = list(
  min_cells = 100L,
  pando_design_args = list(
    peak_to_gene_method = "GREAT",
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

`Pando::prepare_grn_design()` creates a version-2 condition-agnostic dictionary. Exact predictors are deduplicated by TF, measured ATAC feature, and target. The canonical `GREAT` method keeps candidate admission independent of fitted target-expression correlation. A finite `max_edges_per_target` is rejected because candidate order is not an evidence ranking.

## Stage 1 multitask and bootstrap controls

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

Cross-validation remains cell-level and condition-stratified. Bootstrap uses one of two explicit methods:

```text
condition_stratified_sample_cluster_nonparametric
condition_stratified_cell_nonparametric_fallback
```

For sample-cluster bootstrap, each condition contributes its observed number of sample IDs sampled with replacement. Each selected sample contributes all of its cells; total bootstrap cell count may vary. Fallback cell bootstrap retains the original number of cells per condition.

Principal Stage 1 outputs:

```r
step1$grn_result$tf_peak_gene_candidates
step1$grn_result$tf_peak_gene_global
step1$grn_result$tf_peak_gene_condition_all
step1$grn_result$tf_peak_gene_significant
step1$grn_result$condition_target_genes
step1$grn_result$target_model_diagnostics
step1$grn_result$stability_diagnostics
step1$grn_result$celltype_fit_status
step1$grn_result$group_status
step1$grn_result$bootstrap_policy
```

Bootstrap provenance fields include:

```text
bootstrap_method
bootstrap_resampling_unit
bootstrap_sample_col
n_bootstrap_samples_total
min_bootstrap_samples_per_condition
bootstrap_fallback_reason
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

Condition is the RegCompass stratum. Cell type is passed to SuperCell2 as the exact label. The Stage 1 sample column does not alter this geometry.

## Stage 3–5 contracts

A reaction becomes core only when one complete GPR branch is contained in the condition target set. The merged Stage 3 object is a biological reaction catalogue, not a GEM.

```r
layer1_args = list(
  regulatory_alpha = 1,
  gpr_and_method = "min",
  gene_half_saturation = 1
)
```

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

Every condition and metacell under one medium uses the same cached model file, reaction IDs, stoichiometry, and bounds.

## Timing and interpretation

Every stage prints elapsed time and status in the R console. Timing is not stored in returned objects or timing TSV files.

Sample-aware bootstrap measures reproducibility across the observed sample clusters but does not itself constitute donor-level treatment-effect inference. Cross-validation remains cell-level. Metacell-level tests remain descriptive pseudo-observation analyses.

## Documentation

- [Pando and multitask GRN parameter policy](grn-parameter-policy.md)
- [Mathematics and object contracts](multitask-shared-grn.md)
- [Sample-aware bootstrap contract](sample-aware-bootstrap.md)
- [Quick start](tutorial-01-quick-start.md)
- [Stepwise audit](tutorial-02-stepwise-audit.md)
- [Restart and diagnostics](tutorial-03-advanced-restart.md)
- [Targeted reaction remapping](tutorial-04-targeted-reaction-remapping.md)
- [Condition differential analysis](tutorial-05-condition-differential-analysis.md)
- [Stage contracts](stage-interface-contracts.md)
- [Portable execution](portable-execution.md)
