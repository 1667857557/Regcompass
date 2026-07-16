# RegCompassR

RegCompassR runs one supported RNA+ATAC metacell workflow:

```text
condition × sample × cell type
→ SuperCell2 metacells
→ Pando GRN and reaction confidence
→ local meta-module FASTCORE
→ all-strata barrier
→ sample-balanced global gene score and Q95
→ deduplicated global-union GEM
→ metacell-specific penalties
→ directional microCOMPASS scoring
```

## Installation

```r
remotes::install_github("1667857557/SuperCell_Seurat_V4@supercell-2.0")
remotes::install_github("1667857557/Pando_regcompass")
remotes::install_github("1667857557/Regcompass")
```

## Optional setup steps

Prepare or load a RegCompass-compatible GEM before the main workflow. The helper below downloads and converts Human-GEM when you want the packaged default.

```r
library(RegCompassR)

gem <- rc_prepare_human2_gem(version = "2.0.0")
```

Build one shared medium definition for all conditions so downstream scores remain comparable.

```r
medium <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "blood_like"
)
```

## Run RegCompass

Prepare motif PFMs, a genome object and fragment files, then run the main workflow.

```r
library(Pando)
library(BSgenome.Hsapiens.UCSC.hg38)
data(motifs, package = "Pando")

result <- rc_run_regcompass(
  object = object,
  gem = gem,
  outdir = "RegCompass_result",
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  fragment_files = fragment_files,
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type",
  model_mode = "meta_module_gem",
  medium_scenarios = medium,
  metacell_args = list(
    gamma = 150,
    min_cells_per_stratum = 100,
    min_metacell_size = 10,
    min_metacells_per_stratum = 10
  ),
  layer1_args = list(
    local_fastcore = TRUE,
    sample_balance = TRUE,
    expression_batch_correction = "none"
  ),
  pando_args = list(
    min_metacells = 20,
    pando_infer_args = list(
      method = "glm",
      tf_cor = 0.1,
      peak_cor = 0,
      adjust_method = "fdr"
    )
  ),
  layer2_args = list(
    target_direction = "both",
    solver = "highs"
  ),
  upstream_workers = 6,
  layer2_workers = 12,
  parallel_backend = "snow"
)
```

Use one fixed `gamma` for all strict strata. Strata that produce fewer than `pando_args$min_metacells` are recorded as skipped and are not included in downstream global calibration or scoring.

Optional technical-batch correction is configured in `layer1_args` after logCPM merging and before gene scoring:

```r
layer1_args = list(
  expression_batch_correction = "limma",
  technical_batch_cols = "library_batch",
  preserve_design_cols = c("condition", "cell_type")
)
```

Do not use `sample_id` as a removable batch.

Export and downstream statistical summaries can be generated from `result$microcompass` with project-specific code.

## Main result objects

```r
result$layer1$C_rel
result$layer1$reaction_confidence
result$grn_meta_modules$global_reaction_membership
result$microcompass$score
result$microcompass$penalty
result$microcompass$vmax
result$microcompass$feasible
```

See [workflow](docs/workflow.md) and [public functions](docs/functions.md).
