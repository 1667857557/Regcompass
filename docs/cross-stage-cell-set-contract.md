# Stage 1–Stage 2 cell-set contract

Stage 1 applies the fixed `min_cells = 300` filter before normalization and stores the exact analysis-cell IDs in the Seurat object's native cell order:

```r
step1$cell_filter$retained_cells
```

For condition-aware Pando, an undersized condition × cell-type stratum is deleted. A cell type with at least two qualifying conditions uses condition-GRN fitting; a cell type with one qualifying condition uses standard Pando. Cell types without a qualifying stratum remain in `step1$cell_filter$diagnostics` but are not sent to Stage 1 fitting or Stage 2 metacell construction.

Use the Stage 1 result explicitly in a stepwise run:

```r
step1 <- rc_regcompass_step_grn(
  object = A,
  gem = gem,
  outdir = "RegCompass_steps/01_grn",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  condition_col = "dataset",
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

step2 <- rc_regcompass_step_metacells(
  object = A,
  grn = step1,
  outdir = "RegCompass_steps/02_metacells",
  condition_col = "dataset",
  celltype_col = "cell_type"
)
```

With `grn = step1`, Stage 2:

- subsets `A` to exactly the Stage 1 cell-ID set;
- rejects missing or extra cells after subsetting;
- accepts Seurat's native subset order when Seurat does not preserve the requested Stage 1 order;
- rejects conflicting condition, cell-type, RNA-assay, ATAC-assay, or explicit cell-type arguments;
- stores `cell_filter$exact_stage1_match = TRUE` in `step_metacells.rds`;
- records `exact_cell_set`, `order_matches_stage1`, `order_policy`, and `stage1_order_index` under `metacell_object@misc$regcompass_cross_stage_cell_set`.

Omitting `grn` remains supported for backward compatibility. In that mode Stage 2 independently reapplies the same fixed filter, but the result is marked `independent_stage1_filter_reapplication_v1` rather than an exact inherited contract.
