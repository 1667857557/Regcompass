# Public functions in RegCompassR 1.8.8

## Setup and complete runs

- `rc_prepare_gem()`, `rc_prepare_human2_gem()`, `rc_prepare_mouse_gem()`: load bundled pinned GEMs or explicitly rebuild an upstream release.
- `rc_bundled_gem_manifest()`: inspect bundled model provenance.
- `rc_make_medium_scenarios()`: create the shared medium table.
- `rc_run_regcompass()` and `rc_run_regcompass_one_shot()`: execute the six-stage shared-background condition-sub-GRN workflow.
- `rc_parallel_config()`: inspect operating-system backend resolution.

The canonical GRN mode is:

```r
grn_mode = "multitask_shared_backbone"
```

The reproducibility mode for earlier independent Pando fits is:

```r
grn_mode = "legacy_condition_pando"
```

## Inspectable stages

- `rc_regcompass_step_grn()`: build one Pando structural TF–peak–target universe per cell type, then estimate global and condition-specific edge coefficients with condition-balanced elastic net and stability selection.
- `rc_regcompass_step_metacells()`: build condition-level, cell-type-guided SuperCell2 metacells while retaining sample composition.
- `rc_regcompass_step_meta_modules()`: map active condition sub-GRN targets to complete-GPR core reactions and one ordered subsystem/KEGG–Reactome/master-Rhea expansion.
- `rc_regcompass_step_layer1()`: combine RNA support with an ATAC-only projection of the fitted condition sub-GRN and aggregate through GPR rules.
- `rc_regcompass_step_layer2()`: build one medium-specific union GEM, apply one global FASTCORE completion, and run directional COMPASS-like LP scoring.
- `rc_regcompass_step_results()`: assemble GRN provenance, reaction annotations, rankings and condition contrasts.
- `rc_regcompass_step_target_union()`: map additional targets and score them in the exact cached Stage 5 model.

## Stage 1 structural design controls

```r
pando_args = list(
  min_cells = 20L,
  pando_design_args = list(
    peak_to_gene_method = "Signac",
    min_tf_detection = 0,
    min_peak_detection = 0,
    min_target_detection = 0,
    max_edges_per_target = Inf
  )
)
```

`Pando::prepare_grn_design()` creates a condition-agnostic candidate dictionary. Exact design columns are deduplicated by TF, measured ATAC feature and target. `supporting_regions` retains multiple regulatory-region annotations that map to one measured peak.

## Stage 1 multitask controls

```r
multitask_args = list(
  alpha = 0.5,
  global_penalty_factor = 1,
  deviation_penalty_factor = 2,
  lambda_rule = "lambda.1se",
  nfolds = 5L,
  n_stability = 25L,
  stability_fraction = 0.8,
  min_selection_frequency = 0.7,
  min_sign_stability = 0.8,
  min_abs_effect = 0,
  min_cv_rsq = 0,
  candidate_screen_threshold = 0,
  max_edges_per_target = Inf,
  seed = 12345L
)
```

`alpha` must be below one because the symmetric centred condition-deviation representation requires a ridge component for a unique solution. The default candidate screen and edge cap retain the full Pando structural universe.

Principal Stage 1 outputs:

```r
step1$grn_result$tf_peak_gene_candidates
step1$grn_result$tf_peak_gene_global
step1$grn_result$tf_peak_gene_condition_all
step1$grn_result$tf_peak_gene_significant
step1$grn_result$condition_target_genes
step1$grn_result$target_model_diagnostics
step1$grn_result$stability_diagnostics
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

An existing Harmony reduction may replace PCA. Changing cells, assays, reductions, dimensions, seed, gamma or metacell thresholds invalidates the Stage 2 cache.

## Stage 3 condition and merged catalogues

```r
step3$condition_modules$supported_metabolic_genes
step3$condition_modules$core_gene_reaction
step3$condition_modules$reaction_membership
step3$merged_modules$merged_core_reactions
step3$merged_modules$merged_reaction_membership
step3$merged_modules$source_edge_universe_ids
```

A reaction becomes core only when one complete GPR branch is contained in the condition target set. The merged object is a biological reaction catalogue, not a GEM.

## Stage 4 integrated evidence

`layer1_args` accepts:

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

`completion_time_limit` applies only to union-GEM construction. Directional scoring LPs have no time-limit argument. Every condition and metacell under one medium uses the same cached model file, reaction IDs, stoichiometry and bounds.

## Interpretation and plotting

- `rc_build_reaction_annotations()` and `rc_attach_reaction_annotations()`: attach reaction names, formulas, GPRs and database cross-references.
- `rc_test_condition_reactions()`: descriptive same-reaction, same-direction, same-medium comparisons within cell type.
- `rc_report_condition_directions()`: retain forward/reverse targets and report direction-aware support without claiming net flux.
- `rc_select_gene_reactions()`: select scored reactions by GPR gene.
- `rc_plot_condition_reaction()` and `rc_plot_condition_gene_reactions()`: condition plots with reaction and evidence annotations.

The GRN and RNA support are estimated from the same paired multiome dataset. ATAC-only projection and centring reduce direct duplicate weighting but do not create independent validation evidence. Metacell-level tests remain descriptive pseudo-observation analyses.

## Documentation

- [Mathematics and object contracts](multitask-shared-grn.md)
- [Quick start](tutorial-01-quick-start.md)
- [Stepwise audit](tutorial-02-stepwise-audit.md)
- [Restart and diagnostics](tutorial-03-advanced-restart.md)
- [Targeted reaction remapping](tutorial-04-targeted-reaction-remapping.md)
- [Condition differential analysis](tutorial-05-condition-differential-analysis.md)
- [Stage contracts](stage-interface-contracts.md)
- [Portable execution](portable-execution.md)
