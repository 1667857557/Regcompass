# Run modes and stepwise workflow

RegCompass uses one public workflow and resolves the Pando mode automatically.

## Mode selection

```text
condition has >= 2 levels  -> condition_grn
condition has 1 level      -> standard_pando
condition omitted/absent   -> standard_pando
```

### Condition mode

`Pando::infer_condition_grn()` fits the canonical unversioned condition contract
within each broad cell type. The primary penalty uses outer-heldout common-support
TF RNA × peak ATAC projections.

### Standard mode

Original `Pando::infer_grn()` is run within each broad cell type. Standard
TF–peak–gene coefficients are projected to paired cells and aggregated to
metacells. No condition coefficients or condition contrast are calculated.

## One-shot workflow

```r
result <- rc_run_regcompass(
  object = object,
  gem = gem,
  outdir = "result",
  genome = genome,
  condition_col = "condition",  # NULL is allowed
  celltype_col = "cell_type"
)
```

The same entry point handles both modes. `result$analysis_mode` reports the
selected route.

## Stepwise workflow

```r
step1 <- rc_regcompass_step_grn(...)
step2 <- rc_regcompass_step_metacells(...)
step3 <- rc_regcompass_step_meta_modules(...)
step4 <- rc_regcompass_step_layer1(...)
step5 <- rc_regcompass_step_layer2(...)
result <- rc_regcompass_step_results(...)
```

Stage 1 and Stage 2 independently resolve the input design and their
`analysis_mode` values must agree.

## Native SuperCell inputs

Stage 2 passes:

```text
cell.annotation       = broad cell type
cell.split.condition  = condition, or NULL when omitted
gamma                  = requested graining level
```

No combined metadata stratum is created. A metacell is checked after construction
to ensure one broad cell type and, when applicable, one condition.

## Restart boundaries

- Pando mode, motifs, regions, targets, or Pando fitting arguments changed:
  rerun Stage 1 onward.
- reductions, dimensions, gamma, or SuperCell arguments changed: rerun Stage 2
  onward.
- GPR or subsystem annotations changed: rerun Stage 3 onward.
- regulatory support or GPR aggregation changed: rerun Stage 4 onward.
- medium, union GEM, FASTCORE, or LP controls changed: rerun Stage 5 onward.

## Result behavior

Multiple conditions produce condition summaries and pairwise contrasts. A single
effective condition produces reaction rankings and summaries with an empty
condition contrast. Metacell statistics remain within-dataset inference rather
than biological-replicate inference.

See the [public API index](functions.md) and
[stage contracts](stage-interface-contracts.md).
