# Tutorial 2: restartable workflow

Each stage writes a checkpoint to its output directory. Reuse the same input object, GEM, metadata columns, assays, and medium definition when restarting downstream stages.

This tutorial shows only parameters that users commonly need to choose or change. Stable workflow defaults are intentionally omitted.

```r
workers <- 10L
```

`workers` is the workflow-level process cap. Ridge GRNs keep one broad cell type resident at a time and reuse this budget for target-level Pando work.

## 1. Regulatory evidence

Stage 1 selects the Pando route automatically after the `min_cells` filter:

- at least two retained conditions within a broad cell type: condition-comparable common-dictionary ridge;
- one effective condition: standard Pando K=1 ridge;
- `condition_col = NULL`: standard Pando;
- mixed datasets: route independently by broad cell type.

Routine call:

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

Stable defaults include `tf_cor = 0.05`, `peak_cor = 0.05`, BH `padj_threshold = 0.05`, ridge fitting, and `target_rsq_threshold = 0.05`. `target_rsq_threshold` is a RegCompass target-model quality gate on the selected-lambda **final full-data** target R²; it does not redefine Pando edge significance. OOF R² remains diagnostic only.

Only expose changed values. For example:

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
      tf_cor = 0.10,
      peak_cor = 0.10,
      padj_threshold = 0.01
    )
  ),
  target_rsq_threshold = 0.10,
  workers = workers
)
```

## 2. Multimodal metacells

Stage 2 uses the same retained cell/metadata design as Stage 1. For each broad cell type, all conditions jointly build one native RNA+ATAC WNN and Walktrap hierarchy; condition splits membership only **after** clustering.

### 2A. Default path

```r
step2 <- rc_regcompass_step_metacells(
  object = A,
  grn = step1,
  outdir = "run/02_metacells",
  condition_col = step1$params$requested_condition_col,
  celltype_col = step1$params$celltype_col,
  workers = workers
)
```

The default reductions are RNA PCA and ATAC LSI. Only override them when needed:

```r
metacell_args = list(
  rna_reduction = "harmony"
)
```

### 2B. Enforce a minimum final metacell size

When `min_metacell_size > 1`, supply an explicit `min_merge_affinity`. Repair occurs after the condition split and uses only the **original shared WNN**; RegCompass never builds a second condition-specific graph.

```r
step2 <- rc_regcompass_step_metacells(
  object = A,
  grn = step1,
  outdir = "run/02_metacells",
  condition_col = step1$params$requested_condition_col,
  celltype_col = step1$params$celltype_col,
  metacell_args = list(
    min_metacell_size = 10L,
    min_metacells_per_stratum = 3L,
    min_merge_affinity = 0.05,
    unresolved_small_policy = "error"
  ),
  workers = workers
)
```

The numeric affinity threshold is an example, not a hidden biological default. A repair candidate must be in the same condition and broad cell type and satisfy the configured normalized original-WNN affinity. The repaired membership becomes the only canonical membership used by RNA/ATAC aggregation, Pando projection, and downstream scoring.

### 2C. Recount ATAC from fragment files

Fragments are processed only after final repaired membership is fixed.

```r
step2 <- rc_regcompass_step_metacells(
  object = A,
  grn = step1,
  outdir = "run/02_metacells",
  condition_col = step1$params$requested_condition_col,
  celltype_col = step1$params$celltype_col,
  fragment_files = fragment_files,
  workers = workers
)
```

See [metacell-reduction-selection.md](metacell-reduction-selection.md) for the complete geometry, repair, provenance, and fragment contract.

## 3. Reaction meta-modules

```r
step3 <- rc_regcompass_step_meta_modules(
  grn = step1,
  metacells = step2,
  gem = gem,
  outdir = "run/03_meta_modules"
)
```

Only add `meta_module_args$subsystem_table` when the GEM requires an explicit compatible subsystem override.

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

`gpr_and_method` can be changed to `"median"` or `"mean"` when scientifically intended; the default is `"min"`.

## 5. Medium scenarios

Choose the biological medium rather than repeating stable matching controls.

```r
medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "normal_human_plasma",
  species = "human"
)
```

A built-in challenge can be selected similarly:

```r
medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "high_glucose",
  species = "human"
)
```

Built-in challenge concentrations are **concentration metadata**, not automatic uptake-flux constraints. High- and low-concentration presets may therefore resolve to identical LP bounds. To model a flux restriction, supply an explicit custom reaction bound or explicit uptake-fraction assumption. See [medium-presets.md](medium-presets.md).

## 6. Layer 2 structural model and directional scoring

`meta_module_gem` with CORDA2 is the default structural route:

```r
step5 <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "run/05_layer2",
  workers = workers
)
```

Only change CORDA2 controls deliberately. The original controls remain under `model_params$corda2_args`: `MCxNCthresh`, `constraint`, `constrainby`, `om`, and `ci`.

To score only one reaction direction:

```r
step5 <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "run/05_layer2",
  layer2_args = list(target_direction = "forward"),
  workers = workers
)
```

`target_direction` accepts `"both"`, `"forward"`, or `"reverse"`.

FASTCORE remains a supplementary alternative completion route:

```r
step5 <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "run/05_layer2",
  layer2_args = list(
    model_params = list(model_completion = "fastcore")
  ),
  workers = workers
)
```

See [layer2-corda.md](layer2-corda.md) and [layer2-model-builders.md](layer2-model-builders.md) before changing structural controls.

## 7. Result assembly

```r
result <- rc_regcompass_step_results(
  grn = step1,
  metacells = step2,
  meta_modules = step3,
  layer1 = step4,
  layer2 = step5,
  gem = gem,
  outdir = "run/06_results"
)
```

The stage objects are restartable checkpoints. Equations and quantitative definitions are maintained in [mathematical-model.md](mathematical-model.md).
