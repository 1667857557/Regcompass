# Tutorial 2: restartable workflow

Each stage writes a checkpoint. Reuse the same Seurat object, GEM, metadata columns, assays, and medium definition when restarting downstream stages. Only current production parameters are shown here; complete argument definitions are in the corresponding Rd help pages.

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
  pando_args = list(
    min_cells = 500L,
    pando_infer_args = list(
      tf_cor = 0.05,
      peak_cor = 0.05,
      padj_threshold = 0.05,
      rank_action = "mark",
      min_residual_df = 1L,
      reference_condition = "Control"
    )
  ),
  workers = workers
)
```

Use `condition_col = NULL` for a dataset without conditions. Broad cell types with at least two retained conditions use the conditional Pando production route: pooled/global plus condition-local candidate discovery, one frozen exact `(target, TF, peak)` union dictionary, E★ with fixed `z = 0.25`, fusion-component joint-SE inference, BH within each `condition × target`, and an any-condition exact-edge union for RegCompass handoff. Cell types with one effective condition use standard Pando.

For the conditional route, the commonly adjusted Stage 1 parameters are:

- `pando_args$min_cells`: minimum cells required for each retained Stage 1 condition × cell-type stratum.
- `pando_args$pando_infer_args$tf_cor`: TF-target candidate-discovery correlation threshold.
- `pando_args$pando_infer_args$peak_cor`: peak-target candidate-discovery correlation threshold.
- `pando_args$pando_infer_args$padj_threshold`: strict BH threshold after fusion-component joint inference; BH is performed separately within each `condition × target` family.
- `pando_args$pando_infer_args$rank_action`: `"mark"` keeps raw rank-deficient fits and records identifiable/boundary structure; `"error"` requests strict failure.
- `pando_args$pando_infer_args$min_residual_df`: minimum pooled effective residual degrees of freedom.
- `pando_args$pando_infer_args$reference_condition`: predefined experimental reference used to construct the K-condition contrast-tree geometry. It is a design coordinate, not a tuning parameter. The label must be retained in every conditional cell type; if omitted, Pando uses and records the first retained condition. Do not choose it after inspecting GRN results.

The conditional production path does **not** expose `condition_ridge_control`, `condition_e_control`, `cv_folds`, `lambda_rule`, `fusion_ratio`, an alternative `z`, or a sensitivity grid. `z = 0.25` is fixed in the conditional estimator. Target full-data R² is retained as a diagnostic for this route and does not gate the exact-edge handoff.

An exact edge enters RegCompass when every fitted condition has a valid continuous E★ coefficient and at least one condition has `padj < padj_threshold`. Once admitted, that exact edge is kept in every condition with the condition's own continuous `penalty_effect`; a nonsignificant small condition is therefore not converted into an absent edge solely because of lower power.

For a one-condition standard-ridge route, standard Pando retains its separate `ridge_control` API and the existing `target_rsq_threshold` gate. The top-level `target_rsq_threshold` argument is therefore still present for standard Pando; on the conditional E★ route the same R² value is diagnostic only.

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

Stage 2 builds one multimodal WNN graph and one Walktrap hierarchy per broad cell type. Conditions share that graph/hierarchy; final condition-pure metacells are selected by legal cuts of the shared hierarchy. `gamma` is a resolution target, while `min_metacell_size` and `min_metacells_per_stratum` are hard constraints.

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

For conditional GRNs, Layer 1 projects only exact edges admitted by the any-condition BH union, but uses each condition's continuous E★ `penalty_effect`. The existing RegCompass exposure remains

`0.75 × mean(TF × ATAC) + 0.25 × mean(TF) × mean(ATAC)`

within each metacell. This exposure weight `0.25` is independent of the conditional E★ deviation threshold `z = 0.25`.

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

Equations and statistical definitions are maintained in [mathematical-model.md](mathematical-model.md).
