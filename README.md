# RegCompassR

RegCompassR integrates paired single-cell RNA/ATAC regulatory evidence with genome-scale metabolic models for cell-type-resolved reaction analysis.

## Installation

```r
remotes::install_github("1667857557/Pando_regcompass")
remotes::install_github("1667857557/SuperCell_Seurat_V4")
remotes::install_github("1667857557/Regcompass")
```

## Required input

Use a paired-cell Seurat object with RNA and ATAC counts for the same cells, broad cell-type metadata, RNA PCA/Harmony, ATAC LSI, and genome-compatible peak coordinates. `condition_col = NULL` is supported.

## Workflow

```text
paired RNA + ATAC cells
  -> Pando GRN
  -> shared-WNN multimodal metacells
  -> reaction meta-modules
  -> Layer 1 RNA/regulatory evidence
  -> Layer 2 structural model + directional scoring
  -> reaction results and condition comparisons
```

Multi-condition broad cell types use the Pando common-dictionary condition ridge route; one-condition broad cell types use standard Pando. Stage 2 builds one multimodal WNN per broad cell type with conditions jointly present, then splits membership by condition; optional small-metacell repair uses that same original WNN.

## Minimal workflow

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

## Main defaults

| Control | Default |
|---|---|
| Stage 1 `pando_args$min_cells` | `500L` |
| RegCompass target-model `target_rsq_threshold` | `0.05` |
| Stage 2 `min_cells_per_stratum` | `20L` |
| Stage 2 `gamma` / `k.knn` | `30L` / `30L` |
| Structural route | `model_mode = "meta_module_gem"` with CORDA2 |
| Worker cap | `10L` |

When `min_metacell_size > 1L`, also supply `min_merge_affinity`; the final repaired membership is used by all downstream aggregation and scoring.

## Predefined media

`rc_make_medium_scenarios()` supports `normal_human_plasma`, `mouse_plasma`, `high_glucose`, `low_glucose`, `high_lactate`, `low_lactate`, `low_glutamine`, and `custom`. Built-in challenge concentrations are biological concentration metadata and are not automatically converted to uptake flux bounds.

See [medium presets](docs/medium-presets.md) for predefined/custom input formats.

## Documentation

- [Quick start](docs/tutorial-01-quick-start.md)
- [Restartable workflow](docs/tutorial-02-stepwise-audit.md)
- [Function reference](docs/functions.md)
- [Post-analysis](docs/tutorial-04-post-analysis.md)
- [Mathematical specification](docs/mathematical-model.md)
