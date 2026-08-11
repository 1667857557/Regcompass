# Stage 1–Stage 2 cell-set contract

Stage 1 applies `pando_args$min_cells` before normalization and stores the exact analysis-cell IDs in the Seurat object's native cell order. The default is `500L`, but users may supply another positive integer:

```r
step1$cell_filter$retained_cells
```

For condition-aware input, an undersized condition × cell-type stratum is deleted first. Routing is then automatic for each retained broad cell type: at least two retained condition levels use Condition Pando; one retained condition level uses Standard Pando. With `condition_col = NULL`, Standard Pando is used. Mixed datasets may therefore contain both routes in one Stage 1 run.

The stable Pando defaults are `tf_cor = 0.05`, `peak_cor = 0.05`, BH adjustment, and `padj_threshold = 0.05` where the selected route uses that parameter. Known arguments belonging only to the other Pando route are ignored for the incompatible route; genuinely unknown argument names still fail validation.

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
    min_cells = 500L
  )
)

step2 <- rc_regcompass_step_metacells(
  object = A,
  grn = step1,
  outdir = "RegCompass_steps/02_metacells",
  condition_col = "dataset",
  celltype_col = "cell_type",
  metacell_args = list(
    min_cells_per_stratum = 500L
  )
)
```

With `grn = step1`, Stage 2:

- subsets `A` to exactly the Stage 1 cell-ID set;
- rejects missing or extra cells after subsetting;
- accepts Seurat's native subset order when Seurat does not preserve the requested Stage 1 order;
- rejects conflicting condition, cell-type, RNA-assay, ATAC-assay, or explicit cell-type arguments;
- applies the configurable `metacell_args$min_cells_per_stratum` gate, default `500L`, before metacell construction;
- stores `cell_filter$exact_stage1_match = TRUE` in `step_metacells.rds`;
- records `exact_cell_set`, `order_matches_stage1`, `order_policy`, and `stage1_order_index` under `metacell_object@misc$regcompass_cross_stage_cell_set`.

Omitting `grn` remains supported for backward compatibility. In that mode Stage 2 independently reapplies the Stage 1 default `min_cells = 500L` filter, and the result is marked `independent_stage1_filter_reapplication_v1` rather than an exact inherited contract. The Stage 2 `min_cells_per_stratum` default remains separately configurable.