# Tutorial 2: stepwise workflow

Use the stepwise API to save, inspect, and restart individual stages. Equations
are in [Mathematical model](mathematical-model.md).

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

## Stage 1: Pando GRNs

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
      reference_condition = "Control",
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

Inspect:

```r
step1$params$pando_parallel
step1$grn_result$condition_fit_status
step1$grn_result$condition_grn_fits
step1$grn_result$tf_peak_gene_condition
step1$grn_result$tf_peak_gene_condition_effect
step1$grn_result$normalization_policy
```

The absolute condition table is used downstream. The reference-effect table is
for interpretation. See [Pando condition contract](condition-comparable-grn.md).

Mouse runs must supply build-matched regions through
`pando_args$pando_initiate_args$regions`.

## Stage 2: metacells

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

Inspect:

```r
step2$pooled$metacell_meta
step2$pooled$membership
step2$pooled$stratum_status
step2$pooled$cache_contract$analysis_args
```

Set `overwrite = TRUE` after changing cells, assays, reductions, dimensions,
seed, gamma, or thresholds.

## Stage 3: reaction catalogue

```r
step3 <- rc_regcompass_step_meta_modules(
  grn = step1,
  metacells = step2,
  gem = gem,
  outdir = "RegCompass_steps/03_meta_modules"
)
```

Inspect:

```r
step3$condition_modules$supported_metabolic_genes
step3$condition_modules$core_gene_reaction
step3$condition_modules$reaction_membership
step3$merged_modules$merged_core_reactions
step3$merged_modules$merged_reaction_membership
```

Stage 3 requires complete GPR support for core reactions, performs one ordered
annotation expansion, and does not run FASTCORE.

## Stage 4: reaction penalties

```r
step4 <- rc_regcompass_step_layer1(
  grn = step1,
  metacells = step2,
  meta_modules = step3,
  gem = gem,
  outdir = "RegCompass_steps/04_layer1",
  projection_component = "condition",
  comparison_support = "auto",
  regulatory_alpha = 0.5,
  gpr_and_method = "min",
  parallel = TRUE,
  BPPARAM = upstream_bp
)
```

Inspect:

```r
step4$capacity_params
step4$evidence_formula
step4$projection_diagnostics
```

The primary route uses outer-heldout common-support target-gene scores. RNA-only,
condition-full, depth, and alpha routes are sensitivity outputs.

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

Inspect:

```r
step5$model_cache_summary
step5$comparison_table
step5$penalty
step5$vmax
step5$feasible
```

One model is built per medium and reused across conditions and metacells. See
[Medium presets](medium-presets.md).

## Stage 6: results

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

Inspect:

```r
result$reaction_ranking
result$condition_contrast
result$reaction_comparison_by_metacell
result$grn$normalization_policy
```

## Saved stage files

```text
01_grn/step_grn.rds
02_metacells/step_metacells.rds
03_meta_modules/step_meta_modules.rds
04_layer1/step_layer1.rds
05_layer2/step_layer2.rds
06_results/step_results.rds
```

Restart guidance: [Tutorial 3](tutorial-03-advanced-restart.md). Public API:
[functions.md](functions.md).
