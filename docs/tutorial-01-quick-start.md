# Tutorial Level 1: minimal one-shot run

Use this tutorial for a paired-cell RNA+ATAC Seurat object and RegCompassR 1.8.3. See [portable execution](portable-execution.md), [Level 2](tutorial-02-stepwise-audit.md), and [Level 3](tutorial-03-advanced-restart.md) for additional controls.

## Install

```r
install.packages(c("remotes", "highs", "BiocManager"))
BiocManager::install("BiocParallel", ask = FALSE, update = FALSE)
remotes::install_version("SeuratObject", "4.1.4", upgrade = "never")
remotes::install_version("Seurat", "4.4.0", upgrade = "never")
remotes::install_version("Signac", "1.11.0", upgrade = "never")
remotes::install_github("1667857557/SuperCell_Seurat_V4@supercell-2.0", upgrade = "never")
remotes::install_github("1667857557/Pando_regcompass", upgrade = "never")
remotes::install_github("1667857557/Regcompass", upgrade = "never")
```

## Validate input

```r
library(RegCompassR)
library(Pando)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)

data(motifs, package = "Pando")
condition_col <- "dataset"
celltype_col <- "epithelial_or_stem"

stopifnot(
  inherits(A, "Seurat"),
  all(c("RNA", "ATAC") %in% names(A@assays)),
  inherits(A[["ATAC"]], "ChromatinAssay"),
  all(c(condition_col, celltype_col) %in% colnames(A@meta.data)),
  !anyNA(A@meta.data[[condition_col]]),
  !anyNA(A@meta.data[[celltype_col]])
)
```

Use `Pando::motifs` as `pfm`; `motif2tf` is not a motif matrix collection.

## Load the bundled GEM and medium

No model download is required:

```r
gem <- rc_prepare_gem("human")
rc_bundled_gem_manifest()

medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "physiologic",
  species = "human"
)
```

Use `source = "download"` only when intentionally rebuilding or updating the model.

## Run on Linux or Windows

On Linux, set numerical-library threads to one before launching multiple outer workers:

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

```r
result <- rc_run_regcompass_one_shot(
  object = A,
  outdir = "RegCompass_result",
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  fragment_files = FALSE,
  species = "human",
  condition_col = condition_col,
  celltype_col = celltype_col,
  pando_args = list(
    min_cells = 100,
    pando_infer_args = list(
      method = "glm",
      tf_cor = 0.1,
      peak_cor = 0.01,
      adjust_method = "fdr",
      parallel = FALSE
    )
  ),
  metacell_args = list(
    gamma = 75,
    min_cells_per_stratum = 500,
    min_metacell_size = 10
  ),
  layer1_args = list(
    local_fastcore = TRUE,
    local_fastcore_args = list(
      solver = "highs",
      time_limit = 300,
      parallel = TRUE
    )
  ),
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    time_limit = 60
  ),
  upstream_workers = 16L,
  layer2_workers = 12L,
  parallel_backend = "auto",
  progress = TRUE
)
```

`auto` selects SOCK/SnowParam on Windows and MulticoreParam on Linux/macOS. Keep Pando's inner `parallel = FALSE` when outer group parallelism is enabled.

## Confirm completion and timing

```r
stopifnot(
  identical(result$version, "1.8.3"),
  nrow(result$reaction_ranking) > 0,
  file.exists("RegCompass_result/05_layer2/step_layer2.rds"),
  file.exists("RegCompass_result/06_results/regcompass_result.rds"),
  file.exists("RegCompass_result/00_execution_timing.tsv")
)

result$timing$stages
result$timing$total
result$params$parallel_backend_resolved
```

Every stage directory also contains `step_timing.tsv`. Use `progress = FALSE` for quiet batch execution.
