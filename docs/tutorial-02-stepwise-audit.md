# Tutorial 2: restartable workflow

Each stage writes a checkpoint to its output directory. Reuse the same input object, GEM, metadata columns, and medium definition when restarting downstream stages.

```r
workers <- 10L
```

`workers` is the single workflow-level parallel cap. RegCompass selects the platform-specific backend automatically; individual stages use no more workers than their independent task count.

## 1. Regulatory evidence

```r
step1 <- rc_regcompass_step_grn(
  object = A,
  gem = gem,
  outdir = "run/01_grn",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  condition_col = "condition",
  celltype_col = "cell_type",
  pando_args = list(min_cells = 500L),
  workers = workers
)
```

Main parameters are `condition_col`, `celltype_col`, `cell_type`, `rna_assay`, `atac_assay`, `pando_args`, `workers`, and `progress`. Pando inference options, when needed, are supplied through `pando_args$pando_infer_args`; route-specific arguments are dispatched only to the compatible condition or standard Pando path.

## 2. Multimodal metacells

```r
step2 <- rc_regcompass_step_metacells(
  object = A,
  grn = step1,
  outdir = "run/02_metacells",
  condition_col = "condition",
  celltype_col = "cell_type",
  workers = workers
)
```

Main parameters are `fragment_files`, `metacell_args`, `workers`, and `grn`. Supplying `grn = step1` enforces the Stage 1 retained cell set. The public retained-stratum threshold defaults to `500L`; other SuperCell/WNN settings are exposed through `metacell_args` and documented in the Rd help.

## 3. Reaction meta-modules

```r
step3 <- rc_regcompass_step_meta_modules(
  grn = step1,
  metacells = step2,
  gem = gem,
  outdir = "run/03_meta_modules"
)
```

`meta_module_args` contains optional Stage 3 customization, including a custom subsystem table when required.

## 4. Layer 1 reaction evidence

```r
step4 <- rc_regcompass_step_layer1(
  grn = step1,
  metacells = step2,
  meta_modules = step3,
  gem = gem,
  outdir = "run/04_layer1",
  workers = workers
)
```

Public controls are `gpr_and_method`, `gene_half_saturation`, `workers`, and `progress`. The default GPR AND rule is `"min"`.

## 5. Medium scenarios

```r
medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "normal_human_plasma",
  species = "human"
)
```

Built-in scenarios are `normal_human_plasma`, `mouse_plasma`, `high_glucose`, `low_glucose`, `high_lactate`, `low_lactate`, `low_glutamine`, and `custom`. See [medium-presets.md](medium-presets.md) for compositions, provenance, and custom-medium formats.

## 6. Layer 2 structural model and directional scoring

```r
step5 <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "run/05_layer2",
  model_mode = "meta_module_gem",
  workers = workers
)
```

Main parameters are `model_mode`, `layer2_args`, `workers`, and `progress`. `layer2_args` accepts `model_params`, `omega`, `target_direction`, `solver`, and `flux_threshold`. `model_mode = "meta_module_gem"` uses CORDA2 by default; its adjustable `corda2_args` are `MCxNCthresh`, `constraint`, `constrainby`, `om`, and `ci`. Defaults and supplementary routes are documented in the current Rd help and [layer2-corda.md](layer2-corda.md).

## 7. Result assembly

```r
result <- rc_regcompass_step_results(
  grn = step1,
  metacells = step2,
  meta_modules = step3,
  layer1 = step4,
  layer2 = step5,
  gem = gem,
  outdir = "run/06_results",
  species = "human"
)
```

The stage outputs are restartable R objects. For post-analysis functions see [tutorial-04-post-analysis.md](tutorial-04-post-analysis.md) and [functions.md](functions.md). All equations and quantitative definitions are maintained only in [mathematical-model.md](mathematical-model.md).
