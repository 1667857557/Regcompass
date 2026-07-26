# RegCompassR 1.8.8 execution modes

All canonical execution modes follow the same architecture:

```text
one shared Pando structural candidate background per cell type
→ global GRN backbone + condition deviations
→ stability-selected condition sub-GRNs
→ condition metabolic target genes
→ complete-GPR condition cores
→ one ordered biological annotation expansion
→ one shared medium-specific union GEM
→ condition/metacell-specific penalties and directional LP scores
```

## GRN modes

### `grn_mode = "multitask_shared_backbone"`

This is the RegCompassR 1.8.8 default.

- Pando runs once per cell type to build the structural TF–peak–target dictionary.
- All conditions use the same edge universe and edge scaling.
- RegCompass estimates global coefficients and symmetric zero-sum condition deviations.
- Stability selection defines active condition edges.
- `padj` is not used for regularised multitask coefficients.

### `grn_mode = "legacy_condition_pando"`

This mode reproduces independent `condition × cell type` Pando fits. Legacy adjusted-p-value, coefficient and model-R² filters belong only to this mode.

The two GRN modes should be treated as different analyses rather than mixed within one result.

## Level 1: one-shot workflow

Use [Tutorial 1](tutorial-01-quick-start.md) with `rc_run_regcompass_one_shot()`.

- `pando_args` controls structural candidate construction.
- `multitask_args` controls elastic-net fitting, cross-validation and stability selection.
- `upstream_workers` covers Stage 1 and Layer 1.
- `layer2_workers` covers union-GEM construction and LP scoring.

## Level 2: explicit stepwise workflow

Use [Tutorial 2](tutorial-02-stepwise-audit.md) to run:

```r
rc_regcompass_step_grn()
rc_regcompass_step_metacells()
rc_regcompass_step_meta_modules()
rc_regcompass_step_layer1()
rc_regcompass_step_layer2()
rc_regcompass_step_results()
```

Stage 3 maps condition sub-GRN target genes to complete-GPR cores. Stage 3 does not construct a GEM. Stage 5 builds the single medium-specific union GEM used by every condition and metacell.

## Level 3: restart and sensitivity analysis

Use [Tutorial 3](tutorial-03-advanced-restart.md).

- Change motifs, Pando regions, structural detection filters or target genes: rerun Stage 1 onward.
- Change multitask penalties, folds, stability thresholds or candidate screening: rerun Stage 1 onward.
- Change metacell construction: rerun Stage 2 onward.
- Change subsystem annotations or GEM GPR rules: rerun Stage 3 onward.
- Change `regulatory_alpha` or `gpr_and_method`: rerun Stage 4 onward.
- Change medium, FASTCORE or LP settings: rerun Stage 5 onward.

## Level 4: targeted second-pass scoring

Use [Tutorial 4](tutorial-04-targeted-reaction-remapping.md) to remap selected genes or reactions through direct KEGG, Reactome or master-Rhea links and score them in the exact cached Stage 5 union GEM.

## Level 5: condition comparison

Use [Tutorial 5](tutorial-05-condition-differential-analysis.md) to compare scores for the same reaction, direction, medium and cell type between conditions.

## Structural modes

### `model_mode = "meta_module_gem"`

Stage 5 builds one medium-specific union GEM from all condition biological catalogues and performs one global FASTCORE completion.

### `model_mode = "full_gem"`

The complete validated GEM is reused directly. No union reconstruction or FASTCORE completion is performed.

In either structural mode, all conditions within the same analysis use the same stoichiometric matrix and bounds.

See [multitask GRN mathematics and object contracts](multitask-shared-grn.md).
