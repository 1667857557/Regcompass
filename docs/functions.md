# Public functions in RegCompassR 1.8.2

## Setup and complete runs

- `rc_prepare_gem()`, `rc_prepare_human2_gem()`, `rc_prepare_mouse_gem()`: prepare validated species GEMs.
- `rc_make_medium_scenarios()`: create one shared medium table; see [medium presets](medium-presets.md).
- `rc_run_regcompass()` and `rc_run_regcompass_one_shot()`: execute the complete GRN-first workflow.

## Inspectable stages

- `rc_regcompass_step_grn()`: condition-by-cell-type Pando GRNs.
- `rc_regcompass_step_metacells()`: condition-level, cell-type-guided SuperCell2 metacells.
- `rc_regcompass_step_meta_modules()`: complete-GPR cores, subsystem/database expansion, and local FASTCORE completion.
- `rc_regcompass_step_layer1()`: integrated RNA+ATAC reaction expression.
- `rc_regcompass_step_layer2()`: persistent union/full-GEM cache and directional LP scoring.
- `rc_regcompass_step_results()`: rankings, reaction annotations, evidence provenance, and condition contrasts.
- `rc_regcompass_step_target_union()`: use previously scored cores as anchors, expand same-subsystem/KEGG/Reactome/master-Rhea context, and score only annotation-linked reactions that were not global core targets in the original Layer 2 run. FASTCORE-only support and generic union members are excluded.

Stages 3-6 validate workflow parameters, GEM fingerprints, stage classes, and ordered metacell IDs before accepting an upstream object.

## Interpretation and plotting

- `rc_build_reaction_annotations()` and `rc_attach_reaction_annotations()`: reaction names, formulas, GPRs, and evidence classes.
- `rc_test_condition_reactions()`: descriptive same-target comparisons across conditions within cell type.
- `rc_select_gene_reactions()`: select scored reactions by GPR gene.
- `rc_plot_condition_reaction()` and `rc_plot_condition_gene_reactions()`: annotated condition plots.

Sample balancing is not part of the canonical workflow. Metacell-level comparisons are descriptive pseudo-observation analyses rather than automatic biological-replicate inference.

## Tutorials

- [Level 1: quick start](tutorial-01-quick-start.md)
- [Level 2: stepwise audit](tutorial-02-stepwise-audit.md)
- [Level 3: restart and diagnostics](tutorial-03-advanced-restart.md)
- [Expanded non-core target scoring](target-union-scoring.md)
- [Condition statistics and plots](condition-reaction-statistics.md)
