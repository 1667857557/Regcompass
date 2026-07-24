# RegCompassR

RegCompassR 1.8.3 implements a GRN-first RNA+ATAC metabolic workflow:

```text
single-cell RNA normalization
→ cell-type-shared ATAC TF-IDF
→ Pando GRNs for each condition × cell type
→ condition-level, label-guided SuperCell2 metacells
→ complete-GPR core reactions and annotation-expanded meta-modules
→ local FASTCORE completion and one global union GEM
→ RNA+ATAC reaction expression
→ directional COMPASS-like LP scoring
→ optional direct database-linked non-core scoring in the same union GEM
```

## Installation

```r
install.packages(c("remotes", "highs", "BiocManager", "ggplot2"))
BiocManager::install("BiocParallel", ask = FALSE, update = FALSE)
remotes::install_version("SeuratObject", "4.1.4", upgrade = "never")
remotes::install_version("Seurat", "4.4.0", upgrade = "never")
remotes::install_version("Signac", "1.11.0", upgrade = "never")
remotes::install_github("1667857557/SuperCell_Seurat_V4@supercell-2.0", upgrade = "never")
remotes::install_github("1667857557/Pando_regcompass", upgrade = "never")
remotes::install_github("1667857557/Regcompass", upgrade = "never")
```

## Bundled Human-GEM and Mouse-GEM

Human-GEM 2.0.0 and Mouse-GEM 1.8.0 are distributed as validated compressed package assets. Canonical runs therefore require no model download:

```r
human_gem <- rc_prepare_gem("human")
mouse_gem <- rc_prepare_gem("mouse")
rc_bundled_gem_manifest()
```

`source = "bundled"` requires the installed pinned asset and never accesses the network. The download/rebuild path remains available for updates:

```r
updated_gem <- rc_prepare_gem(
  species = "human",
  version = "2.0.0",
  source = "download",
  force_download = TRUE
)
```

The installed manifest records source, release, RDS checksum, size, citation DOI, and CC BY 4.0 attribution. `scripts/build-bundled-gems.R` is the reproducible maintainer build script.

## Input contract

`object` must contain paired RNA and ATAC counts for the same cells, complete condition and cell-type metadata, a Signac `ChromatinAssay`, and peak coordinates matching `genome`. Use `Pando::motifs` as `pfm`; `motif2tf` is not a motif matrix collection.

```r
library(RegCompassR)
library(Pando)
library(BSgenome.Hsapiens.UCSC.hg38)

data(motifs, package = "Pando")

gem <- rc_prepare_gem("human")
medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "physiologic",
  species = "human"
)
```

## Cross-platform parallel execution

`parallel_backend = "auto"` is the recommended setting:

- Windows: `BiocParallel::SnowParam(type = "SOCK")`;
- Linux/macOS: `BiocParallel::MulticoreParam`;
- one worker or unavailable BiocParallel: sequential execution.

```r
rc_parallel_config(workers = 8L, backend = "auto")
```

On Linux, set numerical-library threads to one before starting R when using multiple outer workers:

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

## One-shot run

```r
result <- rc_run_regcompass_one_shot(
  object = A,
  outdir = "RegCompass_result",
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  fragment_files = FALSE,
  species = "human",
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
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
    local_fastcore_args = list(solver = "highs", parallel = TRUE)
  ),
  layer2_args = list(target_direction = "both", solver = "highs"),
  upstream_workers = 16L,
  layer2_workers = 12L,
  parallel_backend = "auto",
  progress = TRUE
)
```

The workflow records the requested and resolved backend, worker counts, and operating system. Use `progress = FALSE` or `options(RegCompassR.progress = FALSE)` for quiet execution.

## Progress and execution time

Every public stage prints a progress indicator and writes `step_timing.tsv` in its output directory. Long BiocParallel loops additionally report task-level progress when RegCompass creates the backend. A complete run additionally writes:

```text
RegCompass_result/00_execution_timing.tsv
```

The final object contains:

```r
result$timing$stages
result$timing$total
result$params$parallel_backend_requested
result$params$parallel_backend_resolved
result$params$upstream_workers
result$params$layer2_workers
```

The timing table reports stage, status, start time, finish time, elapsed seconds, formatted elapsed time, OS type, and R version.

## Score direct database links of selected cores

After a stepwise `meta_module_gem` Layer 2 run, selected cores are mapping anchors only. The function directly identifies non-core reactions sharing a KEGG, Reactome, or master-Rhea ID and runs a second LP only for linked reactions that were not global core targets in the original Layer 2 run.

```r
expanded <- rc_regcompass_step_target_union(
  layer1 = step4,
  meta_modules = step3,
  layer2 = step5,
  gem = gem,
  outdir = "RegCompass_steps/06_expanded_target_scoring",
  core_genes = c("GCLC", "GCLM", "GSS", "GSR", "G6PD", "PGD"),
  gene_match = "complete_gpr",
  progress = TRUE
)
```

Same-subsystem expansion, recursive propagation, FASTCORE-only support, generic union members, metabolite-neighbour reactions, and previously scored global cores are excluded.

## Compare the same reaction across conditions

```r
condition_stats <- rc_test_condition_reactions(
  result,
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
  conditions = c("control_24hr", "JQ1_24hr", "MS177_24hr"),
  cell_types = "stem-cell_like",
  p_adjust_scope = "celltype_contrast_medium"
)

condition_stats$omnibus
condition_stats$pairwise

p <- rc_plot_condition_reaction(
  result,
  reaction_id = "MAR06231",
  cell_type = "stem-cell_like",
  target_direction = "reverse",
  conditions = c("control_24hr", "JQ1_24hr", "MS177_24hr"),
  annotation_p = "p_adj"
)
```

The plot shows one point per metacell and adjusted significance brackets. These are metacell-level, within-dataset comparisons rather than biological-replicate inference. See [Condition-associated reaction statistics](docs/condition-reaction-statistics.md).

## Main outputs

- `01_single_cell_grn/`: GRN status, Pando edges, stage RDS, and timing.
- `02_condition_metacells/`: metacell counts, membership, labels, purity, and timing.
- `03_meta_modules/`: core reactions, expanded modules, FASTCORE diagnostics, and timing.
- `04_layer1/step_layer1.rds`: RNA support, ATAC modifier, GPR diagnostics, reaction expression, and timing.
- `05_layer2/step_layer2.rds`: penalties, relative scores, `vmax`, feasibility, model cache, LP diagnostics, and timing.
- `06_results/regcompass_result.rds`: rankings, annotations, evidence provenance, and timing.
- `00_execution_timing.tsv`: complete stage and total runtime summary.
- optional target-union output: selected core anchors, direct database relation catalog, unique non-core targets, source-model hashes, second-pass LP results, and timing.

## Tutorials

| Level | Use | Tutorial |
|---|---|---|
| 1 | minimal validated one-shot run | [Quick start](docs/tutorial-01-quick-start.md) |
| 2 | stage-by-stage run and audit gates | [Stepwise audit](docs/tutorial-02-stepwise-audit.md) |
| 3 | restart, sensitivity, resources, and failure diagnosis | [Advanced restart](docs/tutorial-03-advanced-restart.md) |

See also [Portable execution, bundled GEMs, progress, and timing](docs/portable-execution.md).
