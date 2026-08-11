# Tutorial 2: restartable workflow

Each stage writes a checkpoint to its output directory. Reuse the same input object, GEM, metadata columns, assays, and medium definition when restarting downstream stages.

This tutorial shows only parameters that users commonly need to choose or change. Stable workflow defaults are intentionally omitted. Use the corresponding Rd help page when you need the complete parameter surface.

```r
workers <- 10L
```

`workers` is the single workflow-level process cap. RegCompass selects the platform-specific backend automatically and shrinks each dispatch to the number of independent tasks available. For ridge GRNs, cell types are processed sequentially and this worker budget is reused inside Pando at the target level; each target worker receives only target-specific multiome data.

## 1. Regulatory evidence

Stage 1 automatically selects the GRN route per retained broad cell type:

- at least two retained conditions: condition-comparable Pando multi-task ridge;
- one effective condition, or `condition_col = NULL`: standard Pando using the K=1 specialization of the same ridge solver;
- mixed datasets: each cell type is routed independently.

`min_cells` is kept visible because it is data dependent. The default correlation thresholds, BH adjustment, ridge CV settings, rank handling, and other stable model defaults are not repeated in the main calls.

### 1A. Condition Pando

Use this when conditions should be compared within each broad cell type.

```r
step1 <- rc_regcompass_step_grn(
  object = A,
  gem = gem,
  outdir = "run/01_grn",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  condition_col = "condition",
  celltype_col = "cell_type",
  pando_args = list(
    min_cells = 500L
  ),
  workers = workers
)
```

Do not set `method` for the condition route. Multi-condition fits use the condition-GRN multi-task ridge engine by design.

Only add `pando_infer_args` when you intentionally want to change a default. For example, a stricter candidate/significance configuration can be supplied as:

```r
pando_args = list(
  min_cells = 500L,
  pando_infer_args = list(
    tf_cor = 0.15,
    peak_cor = 0.10,
    padj_threshold = 0.01
  )
)
```

### 1B. Standard Pando

Use `condition_col = NULL` when the whole analysis should use standard Pando rather than condition-comparable fitting. Ridge is already the default, so `method` should normally be omitted.

```r
step1 <- rc_regcompass_step_grn(
  object = A,
  gem = gem,
  outdir = "run/01_grn",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  condition_col = NULL,
  celltype_col = "cell_type",
  pando_args = list(
    min_cells = 500L
  ),
  workers = workers
)
```

With one effective condition, the shared multi-task ridge backend reduces to K=1: the fusion term is zero, while ridge regularization, target-specific CV, predictor scaling, effective degrees of freedom and ridge-Wald diagnostics remain the same numerical framework.

### 1C. Original standard-Pando GLM compatibility route

Set `method = "glm"` only when the unregularized Gaussian GLM behavior is intentionally required.

```r
step1 <- rc_regcompass_step_grn(
  object = A,
  gem = gem,
  outdir = "run/01_grn",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  condition_col = NULL,
  celltype_col = "cell_type",
  pando_args = list(
    min_cells = 500L,
    pando_infer_args = list(
      method = "glm"
    )
  ),
  workers = workers
)
```

### 1D. Mixed routing

When some broad cell types retain at least two conditions and others retain only one, no model override is required. Condition-capable cell types use the joint multi-task ridge; standard-Pando cell types use its K=1 ridge specialization.

```r
step1 <- rc_regcompass_step_grn(
  object = A,
  gem = gem,
  outdir = "run/01_grn",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  condition_col = "condition",
  celltype_col = "cell_type",
  pando_args = list(
    min_cells = 500L
  ),
  workers = workers
)
```

## 2. Multimodal metacells

Stage 2 should use the same metadata design resolved by Stage 1. The default path aggregates the existing RNA and ATAC count matrices.

### 2A. Existing RNA and ATAC count matrices

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

The default reductions are RNA PCA and ATAC LSI. Only override them when the input object uses different reductions, for example:

```r
metacell_args = list(
  rna_reduction = "harmony"
)
```

### 2B. Recount ATAC from fragment files

Supply `fragment_files` only when metacell ATAC counts should be rebuilt from single-cell fragments.

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

Fragment processing shares the same top-level `workers` cap. Fragment-specific implementation defaults are documented in the Stage 2 Rd page.

## 3. Reaction meta-modules

The default biological expansion rules are workflow contracts and normally require no tuning.

### 3A. Default meta-module construction

```r
step3 <- rc_regcompass_step_meta_modules(
  grn = step1,
  metacells = step2,
  gem = gem,
  outdir = "run/03_meta_modules"
)
```

### 3B. Custom subsystem annotation table

Only add `subsystem_table` when the GEM requires an explicit compatible subsystem annotation override.

```r
step3 <- rc_regcompass_step_meta_modules(
  grn = step1,
  metacells = step2,
  gem = gem,
  outdir = "run/03_meta_modules",
  meta_module_args = list(
    subsystem_table = subsystem_table
  )
)
```

## 4. Layer 1 reaction evidence

The default GPR aggregation and RNA support scaling are normally sufficient.

### 4A. Default Layer 1

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

### 4B. Alternative GPR aggregation

Only expose these arguments when the analysis intentionally changes them.

```r
step4 <- rc_regcompass_step_layer1(
  grn = step1,
  metacells = step2,
  meta_modules = step3,
  gem = gem,
  outdir = "run/04_layer1",
  gpr_and_method = "median",
  workers = workers
)
```

`gpr_and_method` accepts `"min"`, `"median"`, or `"mean"`. `gene_half_saturation` can also be supplied when a non-default RNA support scale is scientifically intended.

## 5. Medium scenarios

Choose the biological medium; do not repeat unchanged exchange and matching defaults in routine calls.

### 5A. Built-in physiological medium

```r
medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "normal_human_plasma",
  species = "human"
)
```

### 5B. Built-in perturbation medium

```r
medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "high_glucose",
  species = "human"
)
```

### 5C. Custom medium

```r
medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "custom",
  species = "human",
  custom_medium = custom_medium
)
```

Built-in scenarios are `normal_human_plasma`, `mouse_plasma`, `high_glucose`, `low_glucose`, `high_lactate`, `low_lactate`, `low_glutamine`, and `custom`. See [medium-presets.md](medium-presets.md) for composition, provenance, and custom-medium formats.

## 6. Layer 2 structural model and directional scoring

`meta_module_gem` with CORDA2 is the default structural route. Stable CORDA2 and scoring defaults are not expanded in the routine call.

### 6A. Default CORDA2 workflow

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

### 6B. Score only one reaction direction

Only set `target_direction` when the analysis should not score both directions.

```r
step5 <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "run/05_layer2",
  layer2_args = list(
    target_direction = "forward"
  ),
  workers = workers
)
```

`target_direction` accepts `"both"`, `"forward"`, or `"reverse"`.

### 6C. Supplementary FASTCORE structural completion

Use this only when intentionally replacing the default CORDA2 completion route.

```r
step5 <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "run/05_layer2",
  layer2_args = list(
    model_params = list(
      model_completion = "fastcore"
    )
  ),
  workers = workers
)
```

CORDA2-specific controls are documented in [layer2-corda.md](layer2-corda.md). FASTCORE-specific structural controls are documented in [layer2-model-builders.md](layer2-model-builders.md).

## 7. Result assembly

Result assembly normally needs no tuning beyond the stage objects, GEM, and output directory.

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

The stage outputs are restartable R objects. Post-analysis functions are listed in [tutorial-04-post-analysis.md](tutorial-04-post-analysis.md) and [functions.md](functions.md). Equations and quantitative definitions are maintained only in [mathematical-model.md](mathematical-model.md).
