# RegCompassR 1.8.4 tutorial index

All execution modes follow the same main workflow:

```text
condition × cell type Pando evidence for Human-GEM genes
→ significantly supported metabolic targets
→ complete-GPR core reactions
→ one ordered subsystem/cross-reference expansion pass
→ integrated RNA+ATAC support with COMPASS GPR-AND aggregation
→ medium-constrained structural model
→ directional LP scoring
```

## Level 1: one-shot workflow

Use [Tutorial 1](tutorial-01-quick-start.md) for a complete run with `rc_run_regcompass_one_shot()`.

- `pfm` is optional; Pando's bundled `motifs` data object is the default.
- `upstream_workers`: Pando inference and Layer 1.
- `layer2_workers`: model completion and LP scoring.
- `meta_module_args`: optional custom subsystem table only.
- `layer1_args`: `regulatory_alpha`, `gpr_and_method`, and gene half-saturation.

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

Stage 3 directly maps significant Pando target genes to complete-GPR cores. It does not perform shared-TF projection or connected-component analysis. Expansion is exactly one ordered pass: core subsystem, direct KEGG/Reactome equivalence, then direct master-Rhea equivalence. The removed interfaces are `expansion_mode`, `max_iterations`, fixed-point recursion, and one-hop reaction expansion.

Stage 4 uses `gpr_and_method = "min"` by default. `"median"` and `"mean"` are available for sensitivity analysis. The former Boltzmann soft-min and `tau` API are removed.

## Level 3: restart and sensitivity analysis

Use [Tutorial 3](tutorial-03-advanced-restart.md).

- Change Pando thresholds, motifs, or regulatory regions: rerun Stage 1 onward.
- Change metacell construction: rerun Stage 2 onward.
- Change subsystem annotations or GEM GPR rules: rerun Stage 3 onward.
- Change `gpr_and_method` or another multiome support setting: rerun Stage 4 onward.
- Change medium, FASTCORE, or LP settings: rerun Stage 5 onward.

## Level 4: targeted second-pass scoring

Use [Tutorial 4](tutorial-04-targeted-reaction-remapping.md) to remap selected genes or core reactions through direct KEGG, Reactome, or master-Rhea links and score the mapped targets in the cached Stage 5 model.

## Level 5: condition comparison

Use [Tutorial 5](tutorial-05-condition-differential-analysis.md) to compare reaction scores between conditions.

## Medium scenarios

Use [Predefined extracellular medium scenarios](medium-presets.md) for physiological, culture-medium, nutrient-sensitivity, technical, and custom options.

## Structural modes

### `model_mode = "meta_module_gem"`

Stage 5 builds the medium-constrained structural model and applies global FASTCORE completion.

### `model_mode = "full_gem"`

The complete validated GEM is reused directly for scoring; no union-model reconstruction or FASTCORE completion is required.

Treat the two modes as separate analyses.
