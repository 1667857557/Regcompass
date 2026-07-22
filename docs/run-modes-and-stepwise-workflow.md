# RegCompassR tutorial index

The tutorials are split by the amount of control and inspection required. Start at the lowest level that matches the analysis task; the biological workflow and canonical defaults remain the same at every level. Examples assume Linux unless stated otherwise.

## Choose a tutorial level

| Level | Use when | What it provides |
|---|---|---|
| [Level 1: minimal one-shot run](tutorial-01-quick-start.md) | the paired-cell Seurat object is ready and a canonical run is needed | installation, minimum input checks, GEM/medium setup, Linux `upstream_workers`/`layer2_workers`, completion check |
| [Level 2: stepwise run with audit gates](tutorial-02-stepwise-audit.md) | every intermediate stage must be inspected | six explicit stages, `MulticoreParam` examples, input/output contracts, files, and continuation gates |
| [Level 3: restart, sensitivity runs, and diagnostics](tutorial-03-advanced-restart.md) | rerunning selected stages or diagnosing structural/resource failures | restart matrix, Linux worker allocation, serial fallbacks, alternative media, full-GEM comparison, solver and Pando diagnostics |

## Canonical workflow shared by all levels

```text
single-cell RNA normalization
→ ATAC TF-IDF shared across conditions within each cell type
→ Pando GRN per condition × cell type (peak_cor = 0.01)
→ cell-type-label-guided SuperCell2 metacells within condition (gamma = 75)
→ post-aggregation membership, purity, and dominant-cell-type audit
→ GRN/metacell coverage validation
→ complete-GPR core reactions
→ subsystem + KEGG/Reactome + master-Rhea expansion
→ local FASTCORE feasibility completion
→ RNA+ATAC reaction expression
→ directional COMPASS-like scoring
```

## Linux parallel map

| Stage | Parallel unit | Main control |
|---|---|---|
| Pando GRN | condition × cell-type group | `upstream_workers` or Stage 1 `BPPARAM` |
| Metacells | no workflow-level BiocParallel unit | metacell arguments |
| Meta-module completion | local FASTCORE per meta-module | `local_fastcore_args$workers` and `$backend` |
| Layer 1 | GPR/reaction capacity | `upstream_workers` or Stage 4 `BPPARAM` |
| Layer 2 | shared-model × metacell LP task | `layer2_workers` or Stage 5 `BPPARAM` |

On Linux, use `parallel_backend = "multicore"` or `BiocParallel::MulticoreParam`. Keep Pando's inner `parallel = FALSE`, and set BLAS/OpenMP thread variables to one before launching R.

## Common input requirements

All levels require:

- a paired-cell Seurat object with RNA counts and an ATAC `ChromatinAssay`;
- RNA, ATAC, and metadata containing the same cell IDs, although order may differ;
- complete condition, cell-type, and selected SuperCell2 label metadata;
- ATAC peak coordinates matching the selected genome;
- a Pando/motifmatchr-compatible PFM/PWM collection such as `Pando::motifs`;
- a supported and validated human or mouse GEM;
- the selected LP solver package.

Do not use the `motif2tf` annotation table as `pfm`. Sample metadata are optional provenance and are not used for balancing, weighting, downsampling, or graph construction.

The canonical label is the column selected by `celltype_col`. One-shot calls can
override it with `metacell_label_col`; Stage 2 calls use `label_col`. The label
guides SuperCell2 before aggregation but does not become an additional hard
stratum, so the post-aggregation composition gates remain required.

## Installation choices

GitHub installation:

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

Local Pando source installation:

```r
install.packages(
  "~/Pando_regcompass.tar.gz",
  repos = NULL,
  type = "source"
)
```

GitHub remote metadata are not required for local Pando installation; the required Pando API is validated at runtime.

## Recommended progression

1. Run Level 1 only after the input contract passes.
2. Use Level 2 for the first real dataset or whenever parameter changes are introduced.
3. Use Level 3 after a valid Level 2 run exists, especially when tuning Linux worker counts or comparing sensitivity runs against an audited baseline.
