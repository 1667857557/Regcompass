# Tutorial Level 1: minimal one-shot run

Use this level when the input object is already a paired-cell RNA+ATAC Seurat object and the goal is to run the canonical workflow with documented defaults.

For stage-by-stage inspection, continue with [Level 2](tutorial-02-stepwise-audit.md). For restart, solver, medium, and resource controls, use [Level 3](tutorial-03-advanced-restart.md).

## 1. Install the required packages

```r
install.packages(c("remotes", "highs", "BiocManager"))
BiocManager::install("BiocParallel", ask = FALSE, update = FALSE)
remotes::install_version("SeuratObject", "4.1.4", upgrade = "never")
remotes::install_version("Seurat", "4.4.0", upgrade = "never")
remotes::install_version("Signac", "1.11.0", upgrade = "never")
remotes::install_github(
  "1667857557/SuperCell_Seurat_V4@supercell-2.0",
  upgrade = "never"
)
remotes::install_github("1667857557/Pando_regcompass", upgrade = "never")
remotes::install_github("1667857557/Regcompass", upgrade = "never")
```

A local Pando source archive is also supported:

```r
install.packages(
  "~/Pando_regcompass.tar.gz",
  repos = NULL,
  type = "source"
)
```

GitHub remote metadata are not required for a local Pando installation. RegCompass validates the Pando API used by the workflow.

## 2. Validate the minimum input contract

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
  !anyNA(A@meta.data[[celltype_col]]),
  all(nzchar(trimws(as.character(A@meta.data[[condition_col]])))),
  all(nzchar(trimws(as.character(A@meta.data[[celltype_col]]))))
)
```

`A` must contain paired RNA and ATAC measurements for the same cell IDs. Assay column order may differ. ATAC peak coordinates and `genome` must use the same genome build. Use the PFM/PWM collection `motifs`; do not use the `motif2tf` annotation table as `pfm`.

## 3. Prepare the GEM and medium

```r
gem <- rc_prepare_gem(
  species = "human",
  version = "2.0.0"
)

medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "high_glucose",
  species = "human"
)

rc_validate_gem(gem)
```

## 4. Run on a Linux multicore system

The following example assumes a Linux host with enough memory for forked workers. `upstream_workers` controls condition × cell-type Pando groups, local FASTCORE completion for individual meta-modules, and Layer 1 GPR-capacity calculations. `layer2_workers` controls the directional LP tasks in Layer 2.

Use fewer Layer 2 workers when each GEM is large because every concurrent solver task needs its own working memory.

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

Keep `pando_infer_args$parallel = FALSE`. RegCompass already distributes the independent condition × cell-type Pando groups across the outer BiocParallel workers; enabling Pando's inner parallelism would create nested workers and CPU oversubscription.

`celltype_col` is automatically passed to SuperCell2 as the pre-aggregation
label. This makes metacell construction label-aware without requiring a second
parameter and reduces cell-type mixing while retaining condition as the only
hard stratum. Purity and the full observed cell-type composition are still
audited after construction.

The analysis order is fixed:

```text
single-cell RNA normalization
→ cell-type-shared ATAC TF-IDF across conditions
→ Pando per condition × cell type [parallel by group]
→ condition-only metacells
→ GRN-derived core reactions and meta-modules
→ local FASTCORE completion [parallel by meta-module]
→ RNA+ATAC reaction expression [parallel by GPR/reaction]
→ directional COMPASS-like scoring [parallel by shared model × metacell]
```

`parallel_backend = "multicore"` explicitly selects `BiocParallel::MulticoreParam` on Linux. Use `"snow"` instead when forked processes are prohibited, such as some containers or managed cluster environments.

## 5. Confirm that the run completed

```r
stopifnot(
  file.exists("RegCompass_result/01_single_cell_grn/pando_group_status.tsv.gz"),
  file.exists("RegCompass_result/02_condition_metacells/metacell_metadata.tsv.gz"),
  file.exists("RegCompass_result/03_meta_modules/core_gene_reaction.tsv.gz"),
  file.exists("RegCompass_result/04_layer1/step_layer1.rds"),
  file.exists("RegCompass_result/05_layer2/step_layer2.rds"),
  file.exists("RegCompass_result/06_results/regcompass_result.rds")
)

head(result$reaction_ranking)
head(result$condition_contrast)
```

Confirm the Stage 3 worker policy:

```r
fastcore_summary <- readRDS(
  "RegCompass_result/03_meta_modules/condition_meta_modules.rds"
)$local_fastcore_summary

unique(fastcore_summary[, c(
  "parallel_task", "parallel_backend", "parallel_workers"
)])
```

Do not interpret the final ranking before confirming that every condition × cell-type Pando group completed and that Layer 2 contains feasible targets. Level 2 shows the required checks at each boundary.
