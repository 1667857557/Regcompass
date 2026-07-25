# RegCompassR 1.8.4 tutorial index

All execution modes use the same architecture:

```text
biological meta-modules
→ merged reaction catalogue
→ Layer 1 multiome support
→ medium-specific union GEM
→ single global FASTCORE
→ directional LP scoring
```

## Level 1: one-shot workflow

Use [Tutorial 1](tutorial-01-quick-start.md) for the complete canonical run through `rc_run_regcompass_one_shot()`.

The one-shot runner manages two worker counts:

- `upstream_workers`: GRN inference and Layer 1;
- `layer2_workers`: union-GEM completion and LP scoring.

Stage 3 constructs the biological catalogue without FASTCORE.

## Level 2: explicit stepwise workflow

Use [Tutorial 2](tutorial-02-stepwise-audit.md) to call:

```r
rc_regcompass_step_grn()
rc_regcompass_step_metacells()
rc_regcompass_step_meta_modules()
rc_regcompass_step_layer1()
rc_regcompass_step_layer2()
rc_regcompass_step_results()
```

Stage 3 returns:

```r
step3$merged_modules$merged_core_reactions
step3$merged_modules$merged_reaction_membership
```

These are catalogue tables, not GEM objects.

## Level 3: restart and sensitivity analysis

Use [Tutorial 3](tutorial-03-advanced-restart.md).

- Change GRN or annotation expansion: rerun Stage 3 onward.
- Change multiome transformation: rerun Stage 4 onward.
- Change medium or global FASTCORE settings: rerun Stage 5 onward.

Global FASTCORE settings are supplied through `layer2_args$model_params`.

## Level 4: targeted second-pass scoring

Use [Tutorial 4](tutorial-04-targeted-reaction-remapping.md) after a completed `meta_module_gem` Stage 5 run.

The second pass reuses the exact final Stage 5 union-GEM files. It verifies the cached file checksum and medium identity, does not rerun FASTCORE, and does not reinterpret the merged Stage 3 catalogue as a GEM.

## Level 5: condition comparison

Use [Tutorial 5](tutorial-05-condition-differential-analysis.md) to compare the same reaction-direction-medium target across conditions.

Within one medium, all conditions share one union GEM. Across different media, structural models may differ because global FASTCORE support may differ.

## Structural modes

### `model_mode = "meta_module_gem"`

- Stage 3: biological catalogue only;
- Stage 5: one medium-specific union GEM;
- one global FASTCORE completion per medium;
- all conditions share the same model within a medium.

### `model_mode = "full_gem"`

- Stage 3 still defines complete-GPR targets;
- the complete validated GEM is reused for scoring;
- no union-GEM reconstruction is performed;
- no FASTCORE reconstruction is required.

Results from these two structural modes are separate analyses and should not be merged into one ranking.
