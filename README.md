# RegCompassR

RegCompassR integrates paired single-cell RNA and ATAC regulatory evidence with cell-type-specific metabolic reaction analysis.

## Installation

```r
remotes::install_github("1667857557/Pando_regcompass")
remotes::install_github("1667857557/SuperCell_Seurat_V4")
remotes::install_github("1667857557/Regcompass")
```

## Required input

Use a paired-cell Seurat object containing:

- RNA and ATAC count assays for the same cells;
- a broad cell-type metadata column;
- RNA PCA (or harmony) and ATAC LSI reductions;
- genome-compatible peak coordinates;
- an optional condition column.

Stage 1 retains groups with at least 300 paired cells. Within each broad cell type, two or more retained conditions use the common-dictionary condition GRN; one retained condition or no condition uses standard Pando automatically.

## Workflow

```text
Paired single-cell RNA + ATAC Seurat object ┐
Species-specific genome-scale model (GEM)  ├──► RegCompass inputs
Shared biological or custom medium          ┘
                                                │
                                                ▼
Stage 1 ── Cell-type Pando GRN routing
          ├── no one condition ──► standard Pando
          └── two or more conditions       ──► common-dictionary condition GRNs
          └──────────────────────────────────► Result: cell-type GRNs, active TF-peak-gene edges, and routing provenance
                                                │
                                                ▼
Stage 2 ── Multimodal metacell construction
          ├── one RNA+ATAC WNN graph per broad cell type
          └── condition split after graph clustering
          └──────────────────────────────────► Result: condition-pure metacells, membership, aggregated RNA/ATAC counts
                                                │
                                                ▼
Stage 3 ── Cell-type reaction meta-modules
          ├── metabolic genes and GPR-linked core reactions
          ├── subsystem / Rhea / network-neighbour expansion
          └──────────────────────────────────► Result: reaction modules and cell-type scoring targets
                                                │
                     Stage 1 GRNs ──────────────┤
                     Stage 2 metacells ─────────┤
                     Stage 3 modules ───────────┘
                                                ▼
Stage 4 ── RNA + regulatory reaction support
          ├── RNA-only gene and reaction support
          └── Pando-regulated multiome support
          └──────────────────────────────────► Result: aligned RNA-only and RNA+ATAC reaction evidence matrices
                                                │
                                                ▼
Stage 5 ── Cell-type × medium structural model
          ├── default       ──► CORDA2 completion
          ├── supplementary ──► FASTCORE completion
          └── supplementary ──► complete full-GEM scoring
          └──────────────────────────────────► Result: shared structural model, bounds, target directions, and directional maximum fluxes
                                                │
                       ┌────────────────────────┴────────────────────────┐
                       ▼                                                 ▼
Primary multiome COMPASS-like penalty                 Matched RNA-only control
                       └────────────────────────┬────────────────────────┘
                                                │
                                                ▼
Stage 6 ── Result assembly and post-analysis
          ├── reaction annotations and evidence classes
          ├── cell-type reaction rankings
          ├── condition contrasts when conditions are available
          └──────────────────────────────────► Result: auditable reaction tables, plots, rankings, and comparisons
```

The primary multiome score and RNA-only control reuse the same completed structural model, bounds, medium, and target directions. Therefore, their difference isolates the regulatory contribution rather than a change in network structure.

## Simple workflow

```r
library(RegCompassR)

result <- rc_run_regcompass_one_shot(
  object = A,
  outdir = "RegCompass_result",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  condition_col = "condition",
  celltype_col = "cell_type",
  pando_args = list(
    min_cells = 300L,
    pando_infer_args = list(
      tf_cor = 0.05,
      peak_cor = 0.05,
      adjust_method = "BH",
      padj_threshold = 0.05,
      rank_action = "mark",
      min_residual_df = 1L
    )
  )
)
```

When no condition metadata are available, omit `condition_col` or set `condition_col = NULL`; the workflow uses standard Pando with one internal background label.

The one-shot workflow prepares the species GEM and default plasma-like medium, builds cell-type-specific GRNs and condition-pure metacells, constructs reaction catalogues, calculates RNA+ATAC-informed reaction penalties, and assembles condition-level results when conditions are available.

Main outputs:

```r
result$grn$cell_type_analysis_mode
result$reaction_catalog
result$reaction_evidence
result$reaction_ranking
result$condition_contrast
result$reaction_comparison_by_metacell
```

## Documentation

1. [Quick run](docs/tutorial-01-quick-start.md)
2. [Stepwise run](docs/tutorial-02-stepwise-audit.md)
3. [Post analysis](docs/tutorial-04-post-analysis.md)
4. [Principles and mathematical formulas](docs/mathematical-model.md)
5. [Function reference](docs/functions.md)
