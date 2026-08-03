# Function index

## Complete workflow

- `rc_run_regcompass_one_shot()`: prepare defaults and run the complete workflow.
- `rc_run_regcompass()`: run the complete workflow with explicit stage arguments.

## Restartable stages

- `rc_regcompass_step_grn()`: filter cells and fit cell-type-specific Pando routes.
- `rc_regcompass_step_metacells()`: construct condition-pure SuperCell metacells.
- `rc_regcompass_step_meta_modules()`: build condition-by-cell-type biological reaction catalogues.
- `rc_regcompass_step_layer1()`: combine RNA and regulatory support and apply GPR rules.
- `rc_regcompass_step_layer2()`: build structural models and score directional reactions.
- `rc_regcompass_step_results()`: assemble annotations, rankings, and contrasts.

## GEM and media

- `rc_prepare_gem()`: load and prepare a supported GEM.
- `rc_validate_gem()`: validate a prepared GEM.
- `rc_make_medium_scenarios()`: construct built-in or custom medium tables.

## Results

- `rc_regcompass_targeted_reactions()`: rerun scoring for a focused reaction list.
- `rc_regcompass_condition_contrast()`: extract condition-level comparisons where available.
- `rc_export_microcompass()`: export Layer 2 matrices and diagnostics.

Use the generated Rd help for complete argument definitions. Mathematical definitions are in [mathematical-model.md](mathematical-model.md).
