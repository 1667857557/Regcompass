# Tutorial Level 1: minimal one-shot run

Use this tutorial for a paired-cell RNA+ATAC Seurat object and RegCompassR 1.8.3. The canonical complete workflow automatically selects the operating-system-specific parallel backend, uses separate upstream and Layer 2 worker counts, and releases every stage worker pool immediately after use.

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

## Configure numerical-library threads

On Linux, set numerical-library threads to one before starting R when multiple outer workers are used:

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

## Run the complete workflow

The defaults are `upstream_workers = 6L`, `layer2_workers = 30L`, and `gamma = 30`. Windows automatically uses SOCK workers; Linux/macOS automatically use multicore workers. Set both worker values to `1L` for a fully serial run.

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
    gamma = 30,
    min_cells_per_stratum = 500,
    min_metacell_size = 10
  ),
  layer1_args = list(
    local_fastcore = TRUE,
    local_fastcore_args = list(
      solver = "highs",
      time_limit = 300
    )
  ),
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    time_limit = 60
  ),
  upstream_workers = 6L,
  layer2_workers = 30L,
  progress = TRUE
)
```

Do not provide `parallel_backend`, `BPPARAM`, or per-FASTCORE worker fields to the canonical complete workflow. Operating-system selection and stage assignment are automatic. Pando's inner `parallel = FALSE` remains required because GRN groups are distributed by the outer upstream worker layer.

Each parallel stage creates, starts, stops, and releases its own pool. Cleanup also runs on errors, and full garbage collection follows pool release. The upstream pool is therefore not retained while Layer 2 is running.

## Confirm completion and timing

```r
stopifnot(
  identical(result$version, "1.8.3"),
  identical(result$params$metacell_gamma, 30L),
  nrow(result$reaction_ranking) > 0,
  file.exists("RegCompass_result/05_layer2/step_layer2.rds"),
  file.exists("RegCompass_result/06_results/regcompass_result.rds"),
  file.exists("RegCompass_result/00_execution_timing.tsv")
)

result$timing$stages
result$timing$total
result$params$parallel_backend_resolved
result$params$upstream_workers
result$params$layer2_workers
result$params$parallel_worker_lifecycle
result$params$parallel_stage_groups
```

Every stage directory also contains `step_timing.tsv`. Use `progress = FALSE` for quiet batch execution.
