# Tutorial 1: quick start

## Input

Use a paired-cell Seurat object containing RNA and ATAC counts for the same cells, broad cell-type metadata, RNA PCA or Harmony, ATAC LSI, and genome-compatible peaks. `condition_col` is optional.

## Prepare the GEM and medium

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

Human workflows default to `normal_human_plasma`; mouse workflows default to `mouse_plasma` when no medium is supplied.

Available biological presets are:

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

Use `scenario = "custom"` with `custom_medium` for reaction-level bounds. Use `scenario = NULL` with `custom_metabolites` for a custom metabolite composition. Exact preset composition, provenance fields such as `background_reference_doi`, `background_validation_reference_doi`, and `challenge_reference_doi`, and custom formats are documented in [medium-presets.md](medium-presets.md).

## Run the workflow

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

### Main arguments

| Argument | Purpose | Default |
|---|---|---|
| `species` | Select human or mouse model defaults | required by one-shot setup |
| `condition_col` | Condition metadata column; use `NULL` when absent | `"condition"` |
| `celltype_col` | Broad cell-type metadata column | `"cell_type"` |
| `rna_assay` / `atac_assay` | Input assay names | `"RNA"` / `"ATAC"` |
| `pando_args$min_cells` | Stage 1 retained-group cell threshold | `500L` |
| `metacell_args$min_cells_per_stratum` | Stage 2 retained-stratum threshold | `500L` |
| `model_mode` | Structural scoring route | `"meta_module_gem"` |
| `workers` | RegCompass-wide worker cap | `10L` |

The two cell-count thresholds are configurable. `model_mode = "meta_module_gem"` uses CORDA2 by default; `model_params$model_completion = "fastcore"` selects the supplementary FASTCORE route, and `model_mode = "full_gem"` selects complete-network scoring. Do not set `model_params$completion_time_limit` for the default CORDA2 route.

### Optional Pando controls

```r
pando_args <- list(
  min_cells = 500L,
  pando_infer_args = list(
    tf_cor = 0.1,
    peak_cor = 0.05,
    adjust_method = "BH",
    padj_threshold = 0.05,
    rank_action = "mark",
    min_residual_df = 1L
  )
)
```

`padj_threshold`, `rank_action`, and `min_residual_df` are condition-GRN controls. Standard Pando receives only arguments supported by its route.

## Main outputs

```r
result$reaction_catalog
result$reaction_evidence
result$reaction_ranking
result$condition_contrast
result$reaction_comparison_by_metacell
```

For restartable stage calls see [tutorial-02-stepwise-audit.md](tutorial-02-stepwise-audit.md). All equations and quantitative definitions are in [mathematical-model.md](mathematical-model.md).
