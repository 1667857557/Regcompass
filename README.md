# RegCompassR

RegCompassR integrates paired single-cell RNA and ATAC regulatory evidence with cell-type-specific metabolic reaction analysis.

## Installation

The package uses the latest default branches of the companion repositories.

```r
remotes::install_github("1667857557/Pando_regcompass")
remotes::install_github("1667857557/SuperCell_Seurat_V4")
remotes::install_github("1667857557/Regcompass")
```

## Required input

The input is a paired-cell Seurat object containing:

- RNA and ATAC count assays for the same cells;
- a broad cell-type metadata column;
- RNA PCA and ATAC LSI reductions;
- genome-compatible peak coordinates;
- an optional condition column.

Stage 1 retains condition-by-cell-type strata . Within each cell type, two or more retained conditions use the common-dictionary condition GRN; one retained condition uses standard Pando automatically.

## One-shot workflow

```r
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
      tf_cor = 0.1,
      peak_cor = 0,
      adjust_method = "BH",
      padj_threshold = 0.05,
      rank_action = "mark",
      min_residual_df = 1L
    )
  ),
  metacell_args = list(
    rna_reduction = "pca",
    rna_dims = 1:30,
    atac_reduction = "lsi",
    atac_dims = 2:30,
    gamma = 30L,
    min_cells_per_stratum = 300L,
    min_metacell_size = 10L,
    min_metacells_per_stratum = 2L
  )
)
```

## Documentation

- [Quick start](docs/tutorial-01-quick-start.md)
- [Restartable stages](docs/tutorial-02-stepwise-audit.md)
- [Targeted reaction remapping](docs/tutorial-04-targeted-reaction-remapping.md)
- [Condition-level results](docs/tutorial-05-condition-differential-analysis.md)
- [Mathematical specification](docs/mathematical-model.md)
- [Function index](docs/functions.md)
