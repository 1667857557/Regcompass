# Tutorial 2: restartable workflow

Each stage writes a checkpoint. Reuse the same Seurat object, GEM, metadata columns, assays, and medium definition when restarting downstream stages. Only commonly adjusted parameters are shown here; complete argument definitions are in the corresponding Rd help pages.

```r
workers <- 10L
```

## 1. Regulatory GRN

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

Use `condition_col = NULL` for a dataset without conditions. RegCompass automatically uses condition-comparable Pando ridge when a retained cell type has at least two conditions and standard Pando ridge otherwise.

Only specify `pando_infer_args` when changing a validated default, for example:

```r
pando_args = list(
  min_cells = 500L,
  pando_infer_args = list(
    tf_cor = 0.10,
    peak_cor = 0.10,
    padj_threshold = 0.01
  )
)
```

## 2. Multimodal metacells

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

Stage 2 builds one multimodal WNN graph and one Walktrap hierarchy per broad cell type. Conditions share that graph/hierarchy; final condition-pure metacells are selected from the shared hierarchy. `gamma` is a resolution target, while `min_metacell_size` and `min_metacells_per_stratum` are hard constraints.

Common overrides:

```r
metacell_args = list(
  rna_reduction = "harmony",
  gamma = 30L,
  min_metacell_size = 5L,
  min_metacells_per_stratum = 2L
)
```

Supply `fragment_files` to `rc_regcompass_step_metacells()` only when metacell ATAC counts should be rebuilt from fragments.

## 3. Reaction meta-modules

```r
step3 <- rc_regcompass_step_meta_modules(
  grn = step1,
  metacells = step2,
  gem = gem,
  outdir = "run/03_meta_modules"
)
```

Use `meta_module_args = list(subsystem_table = subsystem_table)` only for an intentional compatible subsystem override.

## 4. Layer 1 regulatory reaction evidence

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

`gpr_and_method` accepts `"min"`, `"median"`, or `"mean"`. Quantitative RNA is computed from single-cell linear CPM and averaged equally within the exact final SuperCell membership.

## 5. Medium

Built-in scenarios:

- `normal_human_plasma`
- `mouse_plasma`
- `high_glucose`
- `low_glucose`
- `high_lactate`
- `low_lactate`
- `low_glutamine`
- `custom`

```r
medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "normal_human_plasma",
  species = "human"
)
```

For a custom medium:

```r
medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "custom",
  species = "human",
  custom_medium = custom_medium
)
```

See [medium-presets.md](medium-presets.md) for predefined-medium composition/provenance and custom table requirements.

## 6. Layer 2

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

`meta_module_gem` with CORDA2 is the default structural route. Change `layer2_args` only when a different target direction, completion route, or CORDA2 control is intentionally required.

Example:

```r
layer2_args = list(target_direction = "forward")
```

## 7. Results

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

Equations and statistical definitions are maintained in [mathematical-model.md](mathematical-model.md), not duplicated in the tutorial.
