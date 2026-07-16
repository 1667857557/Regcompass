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

## 1. Prepare Human-GEM

```r
library(RegCompassR)

gem <- rc_prepare_human2_gem(version = "2.0.0")
```

## 2. Define a shared medium

```r
medium <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "blood_like"
)
```

The integrated workflow requires the same medium constraints for all conditions.

## 3. Run RegCompass

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
    adaptive_gamma = TRUE,
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

Optional technical-batch correction is configured in `layer1_args`:

```r
layer1_args = list(
  expression_batch_correction = "limma",
  technical_batch_cols = "library_batch",
  preserve_design_cols = c("condition", "cell_type")
)
```

Do not use `sample_id` as a removable batch.

## 4. Test condition differences

```r
differential <- rc_test_microcompass_differential(
  result = result$microcompass,
  formula = score ~ condition,
  method = "limma_continuous"
)
```

Meta cells are aggregated to biological samples before testing.

## 5. Export matrices and diagnostics

```r
rc_export_microcompass(
  result = result$microcompass,
  outdir = "RegCompass_result/export"
)
```

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
