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

When the top-level workflow creates the medium automatically, human models use `normal_human_plasma` and mouse models use `mouse_plasma`.

Built-in scenarios are `normal_human_plasma`, `mouse_plasma`, `high_glucose`, `low_glucose`, `high_lactate`, `low_lactate`, `low_glutamine`, and `custom`. See [medium-presets.md](medium-presets.md) for compositions, provenance, and custom-medium formats.

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
| `species` | Select human or mouse model defaults | required |
| `condition_col` | Condition metadata column; use `NULL` when absent | `"condition"` |
| `celltype_col` | Broad cell-type metadata column | `"cell_type"` |
| `rna_assay` / `atac_assay` | Input assay names | `"RNA"` / `"ATAC"` |
| `pando_args` | Stage 1 options; `min_cells` is the retained-group threshold | `min_cells = 500L` |
| `metacell_args` | Stage 2 SuperCell/WNN options | public retained-stratum threshold `500L` |
| `model_mode` | Layer 2 structural route | `"meta_module_gem"` |
| `layer2_args` | Layer 2 scoring/reconstruction options | `list()` |
| `workers` | RegCompass-wide worker cap | `10L` |

`model_mode = "meta_module_gem"` uses CORDA2 by default. Advanced stage-specific options are documented by the corresponding Rd help pages rather than duplicated in this tutorial.

## Main outputs

```r
result$reaction_catalog
result$reaction_evidence
result$reaction_ranking
result$condition_contrast
result$reaction_comparison_by_metacell
```

For restartable calls see [tutorial-02-stepwise-audit.md](tutorial-02-stepwise-audit.md). All equations and quantitative definitions are maintained only in [mathematical-model.md](mathematical-model.md).
