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
  pando_args = list(
    min_cells = 500L,
    pando_infer_args = list(
      tf_cor = 0.1,
      peak_cor = 0.05,
      adjust_method = "BH",
      padj_threshold = 0.05,
      rank_action = "mark",
      min_residual_df = 1L
    )
  ),
  workers = workers
)
```

Key arguments: `condition_col`, `celltype_col`, `cell_type`, `rna_assay`, `atac_assay`, `pando_args`, `workers`, and `progress`. `pando_args$min_cells` defaults to `500L`. Condition-only Pando arguments are routed only to the condition-GRN path.

## 2. Multimodal metacells

```r
step2 <- rc_regcompass_step_metacells(
  object = A,
  grn = step1,
  outdir = "run/02_metacells",
  condition_col = "condition",
  celltype_col = "cell_type",
  metacell_args = list(
    rna_reduction = "pca",
    rna_dims = 1:30,
    atac_reduction = "lsi",
    atac_dims = 2:30,
    gamma = 30L,
    min_cells_per_stratum = 500L,
    min_metacell_size = 10L,
    min_metacells_per_stratum = 2L
  ),
  workers = workers
)
```

Key arguments: `fragment_files`, `metacell_args`, `workers`, and `grn`. Supplying `grn = step1` enforces the Stage 1 cell set. `metacell_args$min_cells_per_stratum` defaults to `500L`. Fragment aggregation, when requested, uses the same top-level `workers` cap; `metacell_args$fragment_args$workers` is not a public control.

## 3. Reaction meta-modules

```r
step3 <- rc_regcompass_step_meta_modules(
  grn = step1,
  metacells = step2,
  gem = gem,
  outdir = "run/03_meta_modules"
)
```

`meta_module_args` is an optional list for Stage 3 customization, including a custom subsystem table when required.

## 4. Layer 1 reaction evidence

```r
step4 <- rc_regcompass_step_layer1(
  grn = step1,
  metacells = step2,
  meta_modules = step3,
  gem = gem,
  outdir = "run/04_layer1",
  gpr_and_method = "min",
  gene_half_saturation = 1,
  workers = workers
)
```

`gpr_and_method` accepts `"min"`, `"median"`, or `"mean"`; the default is `"min"`. `gene_half_saturation` controls the bounded structural-support path. Mathematical definitions are not duplicated here; see [mathematical-model.md](mathematical-model.md).

## 5. Medium scenarios

```r
medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "normal_human_plasma",
  species = "human"
)
```

Built-in scenarios are:

```text
normal_human_plasma
mouse_plasma
high_glucose
low_glucose
high_lactate
low_lactate
low_glutamine
custom
```

Multiple compatible scenarios may be supplied as a character vector. Use `scenario = "custom"` with `custom_medium` for reaction-level bounds, or `scenario = NULL` with `custom_metabolites` for a custom metabolite composition. Provenance fields including `background_reference_doi`, `background_validation_reference_doi`, and `challenge_reference_doi` are retained in the returned table. See [medium-presets.md](medium-presets.md).

## 6. Layer 2 structural model and directional scoring

```r
step5 <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "run/05_layer2",
  model_mode = "meta_module_gem",
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    model_params = list(
      strict = TRUE,
      corda2_args = list(
        MCxNCthresh = 2,
        constraint = 1,
        constrainby = "val",
        om = 1e4,
        ci = 0.01
      )
    )
  ),
  workers = workers
)
```

`layer2_args` accepts `model_params`, `omega`, `target_direction`, `solver`, and `flux_threshold`. `model_mode = "meta_module_gem"` uses CORDA2 by default. Set `model_params$model_completion = "fastcore"` for supplementary FASTCORE completion or `model_mode = "full_gem"` for complete-network scoring. `model_params$completion_time_limit` is not accepted for CORDA2.

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

The stage outputs are restartable R objects. For post-analysis functions see [tutorial-04-post-analysis.md](tutorial-04-post-analysis.md) and [functions.md](functions.md). All equations and algorithmic definitions are maintained in [mathematical-model.md](mathematical-model.md).
