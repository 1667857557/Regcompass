# Tutorial 2: stepwise workflow

Use the stepwise API to save, inspect and restart the canonical stages.
Equations are in [Tutorial 3](tutorial-03-mathematical-model.md). Public API:
[functions.md](functions.md).

## Parallel backends

```r
library(BiocParallel)

upstream_bp <- if (.Platform$OS.type == "windows") {
  SnowParam(workers = 6L, type = "SOCK", progressbar = TRUE)
} else {
  MulticoreParam(workers = 6L, progressbar = TRUE)
}

layer2_bp <- if (.Platform$OS.type == "windows") {
  SnowParam(workers = 30L, type = "SOCK", progressbar = TRUE)
} else {
  MulticoreParam(workers = 30L, progressbar = TRUE)
}
```

`BPPARAM = TRUE` is invalid. Do not set `parallel` inside
`pando_infer_args`.

## Stage 1: condition-aware or standard Pando

```r
step1 <- rc_regcompass_step_grn(
  object = A,
  gem = gem,
  outdir = "RegCompass_steps/01_grn",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  condition_col = "Group",
  celltype_col = "cell_type",
  pando_args = list(
    min_cells = 100L,
    min_abs_estimate = 0,
    min_model_rsq = 0.1,
    pando_infer_args = list(
      candidate_screen = "motif_domain",
      condition_mix = 0.5,
      condition_weight = "equal",
      outer_nfolds = 5L,
      inner_nfolds = 5L,
      lambda_selection = "lambda.1se",
      scale = TRUE
    )
  ),
  parallel = TRUE,
  BPPARAM = upstream_bp
)
```

```r
step1$params$analysis_mode
step1$grn_result$condition_grn_fits
step1$grn_result$tf_peak_gene_condition
```

With fewer than two condition levels, Stage 1 uses `standard_pando` and No
condition coefficients are calculated.

## Stage 2: cell-type graphs and condition-pure metacells

```r
step2 <- rc_regcompass_step_metacells(
  object = A,
  outdir = "RegCompass_steps/02_metacells",
  condition_col = "Group",
  celltype_col = "cell_type",
  fragment_files = FALSE,
  metacell_args = list(
    rna_reduction = "pca",
    rna_dims = 1:30,
    atac_reduction = "lsi",
    atac_dims = 2:30,
    gamma = 30L,
    seed = 12345L,
    min_cells_per_stratum = 500L,
    min_metacell_size = 10L,
    min_metacells_per_stratum = 2L,
    overwrite = FALSE
  )
)
```

```r
step2$pooled$metacell_meta
step2$pooled$membership
step2$pooled$input_design
```

Each cell type has one independent graph. All conditions are joint in that
graph, and final metacells are condition-pure. Set `overwrite = TRUE` after
changing cells, reductions, dimensions, gamma, seed or thresholds.

## Stage 3: reaction catalogue

```r
step3 <- rc_regcompass_step_meta_modules(
  grn = step1,
  metacells = step2,
  gem = gem,
  outdir = "RegCompass_steps/03_meta_modules"
)
```

```r
step3$condition_modules$core_gene_reaction
step3$merged_modules$merged_core_reactions
step3$merged_modules$merged_reaction_membership
```

## Stage 4: condition-full reaction support

```r
step4 <- rc_regcompass_step_layer1(
  grn = step1,
  metacells = step2,
  meta_modules = step3,
  gem = gem,
  outdir = "RegCompass_steps/04_layer1",
  projection_component = "condition",
  comparison_support = "auto",
  regulatory_alpha = 1,
  gpr_and_method = "min",
  parallel = TRUE,
  BPPARAM = upstream_bp
)
```

```r
step4$gene_projection_condition_full_oof
step4$gene_projection_common_oof
step4$gene_projection_condition_unique_oof
step4$reaction_expression_condition_full_oof
step4$reaction_expression_common_oof
step4$projection_provenance
```

Condition-full OOF is primary. Common support is the jointly estimable
component. Each non-estimable edge side contributes zero.

## Stage 5: shared model and LP scoring

```r
step5 <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "RegCompass_steps/05_layer2",
  model_mode = "meta_module_gem",
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    model_params = list(
      completion_time_limit = 600,
      fastcore_epsilon = 1e-4,
      max_support_reactions = 2000,
      strict = TRUE
    )
  ),
  parallel = TRUE,
  BPPARAM = layer2_bp
)
```

```r
step5$penalty_condition_full_oof
step5$penalty_common_oof
step5$penalty_condition_unique_increment
step5$penalty_rna_only
step5$model_cache_summary
```

All routes use the exact same medium-specific model, bounds, target directions
and `vmax`. The five retired guardrails are absent from the result schema.

## Stage 6: final result

```r
result <- rc_regcompass_step_results(
  grn = step1,
  metacells = step2,
  meta_modules = step3,
  layer1 = step4,
  layer2 = step5,
  gem = gem,
  outdir = "RegCompass_steps/06_results",
  species = "human"
)
```

```r
result$reaction_ranking
result$condition_contrast
result$common_support_component_summary
result$condition_unique_penalty_increment_summary
```

Saved stages:

```text
01_grn/step_grn.rds
02_metacells/step_metacells.rds
03_meta_modules/step_meta_modules.rds
04_layer1/step_layer1.rds
05_layer2/step_layer2.rds
06_results/regcompass_result.rds
```

Use [Tutorial 4](tutorial-04-condition-differential-analysis.md) for condition
statistics.
