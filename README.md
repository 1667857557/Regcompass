# RegCompassR

RegCompassR integrates paired single-cell RNA and ATAC regulatory evidence with genome-scale metabolic models and returns cell-type-resolved reaction scores. The scores are model-derived reaction support/penalty measures, not direct flux measurements.

## Installation

```r
remotes::install_github("1667857557/Pando_regcompass")
remotes::install_github("1667857557/SuperCell_Seurat_V4")
remotes::install_github("1667857557/Regcompass")
```

## Required input

Use a paired-cell Seurat object with RNA and ATAC count assays for the same cells, broad cell-type metadata, RNA PCA or Harmony, ATAC LSI, and genome-compatible peak coordinates. A condition column is optional.

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

For condition-aware analysis, cell types with at least two retained conditions use the condition-GRN route; cell types with one retained condition use standard Pando. `condition_col = NULL` is valid for analyses without condition metadata.

## Main defaults

- Human GEM: Human-GEM `2.0.0`; mouse GEM: Mouse-GEM `1.8.0`.
- Human default medium: `normal_human_plasma`; mouse default medium: `mouse_plasma`.
- Stage 1: `pando_args$min_cells = 500L`.
- Pando candidate defaults: `tf_cor = 0.1`, `peak_cor = 0.05`, `adjust_method = "BH"`; the condition route additionally uses `padj_threshold = 0.05`, `rank_action = "mark"`, and `min_residual_df = 1L`.
- Stage 2: `metacell_args$min_cells_per_stratum = 500L`.
- Structural mode: `model_mode = "meta_module_gem"`, with CORDA2 as the default completion method. `model_params$completion_time_limit` is not a CORDA2 control; it is reserved for supplementary non-CORDA2 completion such as FASTCORE.
- Parallelism: one top-level `workers` cap, default `10L`.

The cell-count thresholds are defaults rather than fixed requirements and can be overridden explicitly.

## Medium presets

`rc_make_medium_scenarios()` supports:

`normal_human_plasma`, `mouse_plasma`, `high_glucose`, `low_glucose`, `high_lactate`, `low_lactate`, `low_glutamine`, and `custom`.

See [`docs/medium-presets.md`](docs/medium-presets.md) for compositions, provenance, and custom-medium input formats.

## Documentation

- [Quick start](docs/tutorial-01-quick-start.md)
- [Restartable stepwise workflow](docs/tutorial-02-stepwise-audit.md)
- [Post analysis](docs/tutorial-04-post-analysis.md)
- [Function reference](docs/functions.md)
- [Mathematical specification](docs/mathematical-model.md)

Algorithmic equations and quantitative definitions are maintained only in `docs/mathematical-model.md`; tutorials and Rd files are interface documentation.
