# Tutorial 1: quick start

## Input

Use a paired-cell Seurat object with RNA and ATAC counts for the same cells, broad cell-type metadata, an RNA reduction (`pca` by default), an ATAC LSI reduction, and genome-compatible peaks. `condition_col` may be `NULL`.

## GEM and predefined medium

```r
library(RegCompassR)

gem <- rc_prepare_gem(
  species = "human",
  version = "2.0.0",
  source = "bundled"
)

medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "normal_human_plasma",
  species = "human"
)
```

Built-in scenarios are `normal_human_plasma`, `mouse_plasma`, `high_glucose`, `low_glucose`, `high_lactate`, `low_lactate`, `low_glutamine`, and `custom`. If omitted, the workflow uses `normal_human_plasma` for human and `mouse_plasma` for mouse. Challenge concentrations are metadata unless an explicit flux/bound assumption is supplied.

## Run

```r
result <- rc_run_regcompass_one_shot(
  object = A,
  outdir = "RegCompass_result",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  gem = gem,
  medium_scenarios = medium_scenarios,
  condition_col = "condition",
  celltype_col = "cell_type",
  workers = 10L
)
```

## Main arguments

| Argument | Purpose | Default |
|---|---|---|
| `condition_col` | Condition metadata; `NULL` selects standard Pando | `"condition"` |
| `celltype_col` | Broad cell-type metadata | `"cell_type"` |
| `rna_assay` / `atac_assay` | Input assays | `"RNA"` / `"ATAC"` |
| `pando_args` | Stage 1 Pando options | `min_cells = 500L` |
| `target_rsq_threshold` | RegCompass target-model final-fit R² gate | `0.05` |
| `metacell_args` | Stage 2 WNN/metacell options | `min_cells_per_stratum = 20L` |
| `model_mode` | Layer 2 route | `"meta_module_gem"` |
| `layer2_args` | Layer 2 options | `list()` |
| `workers` | Workflow worker cap | `10L` |

`model_mode = "meta_module_gem"` uses CORDA2 by default. Use the Rd pages for the complete parameter surface.

## Main outputs

```r
result$reaction_catalog
result$reaction_evidence
result$reaction_ranking
result$condition_contrast
result$reaction_comparison_by_metacell
```

See [tutorial-02-stepwise-audit.md](tutorial-02-stepwise-audit.md) for restartable stage calls and [medium-presets.md](medium-presets.md) for predefined/custom media.
