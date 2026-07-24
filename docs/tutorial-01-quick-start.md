# Tutorial Level 1: minimal one-shot run

Use this tutorial for a paired-cell RNA+ATAC Seurat object and the canonical RegCompassR 1.8.2 workflow. Use [Level 2](tutorial-02-stepwise-audit.md) for stage inspection and [Level 3](tutorial-03-advanced-restart.md) for restart and diagnostics.

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

A local Pando source archive is valid:

```r
install.packages("~/Pando_regcompass.tar.gz", repos = NULL, type = "source")
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

Use `Pando::motifs` as `pfm`; `motif2tf` is not a motif matrix collection. RNA and ATAC must represent the same cell IDs, although assay column order may differ.

## Prepare GEM and medium

```r
gem <- rc_prepare_gem(species = "human", version = "2.0.0")
medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "physiologic",
  species = "human"
)
```

See [medium presets](medium-presets.md) before replacing the physiological baseline. Medium constraints never create a reaction direction absent from the source GEM.

## Run on Linux

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

```r
upstream_workers <- 16L
layer2_workers <- 12L

result <- rc_run_regcompass_one_shot(
  object = A,
  outdir = "RegCompass_result",
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  fragment_files = FALSE,
  species = "human",
  gem = gem,
  medium_scenarios = medium_scenarios,
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
  upstream_workers = upstream_workers,
  layer2_workers = layer2_workers,
  parallel_backend = "multicore"
)
```

Keep Pando's inner `parallel = FALSE`. `celltype_col` is automatically passed to SuperCell2 as the construction label, while condition remains the only hard metacell stratum.

## Confirm completion

```r
stopifnot(
  identical(result$version, "1.8.2"),
  identical(result$schema_version, "regcompass_grn_first_v2"),
  nrow(result$reaction_ranking) > 0,
  nrow(result$reaction_catalog) > 0,
  nrow(result$reaction_evidence) > 0,
  file.exists("RegCompass_result/05_layer2/step_layer2.rds"),
  file.exists("RegCompass_result/06_results/regcompass_result.rds")
)
```

The one-shot output includes stage classes, GEM provenance, model-cache diagnostics, LP penalties, and reaction annotations. Expanded target scoring requires the classed stepwise Stage 3-5 objects described in Level 2.
