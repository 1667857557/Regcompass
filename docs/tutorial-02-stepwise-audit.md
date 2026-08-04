# Tutorial 2: restartable workflow

Each stage writes an RDS checkpoint to its output directory.

## Stage 1: regulatory evidence

```r
library(BiocParallel)

upstream_bp <- if (.Platform$OS.type == "windows") {
  SnowParam(workers = 6L, type = "SOCK", progressbar = TRUE)
} else {
  MulticoreParam(workers = 6L, progressbar = TRUE)
}

step1 <- rc_regcompass_step_grn(
  object = A,
  gem = gem,
  outdir = "run/01_grn",
  genome = BSgenome.Hsapiens.UCSC.hg38,
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
  parallel = TRUE,
  BPPARAM = upstream_bp
)
```

Stage 1 resolves the route independently for each retained cell type:

- at least two retained conditions: common-dictionary condition GRN;
- one retained condition: standard Pando;
- no retained stratum: excluded.

The same `pando_infer_args` list can be used for all routes. Condition-only arguments (`padj_threshold`, `rank_action`, `min_residual_df`, and layer controls) are disabled before standard Pando is called. Standard-model controls (`method`, `alpha`, `scale`, and related model arguments) are disabled for the fixed common-dictionary condition model. Unknown argument names still fail before model fitting.

`tf_cor = 0.1` controls Pando candidate discovery. It must not be interpreted as the final RegCompass penalty-entry threshold. In both the condition-GRN and standard-Pando routes, an edge enters the regulatory penalty only when:

```r
estimable == TRUE
padj < 0.05
abs(corr) >= 0.05
abs(estimate) >= 0.05
```

The absolute correlation and effect-size cutoffs are inclusive. A coefficient with `corr = -0.05` and `estimate = -0.05` passes those two gates; values with absolute magnitude below `0.05` do not. Edges that fail a gate are retained for audit in the complete table but use `penalty_effect = 0`.

For common-dictionary condition fits, `corr_source` records whether `corr` came directly from the condition coefficient table or from the frozen dictionary's `max_abs_tf_target_cor`. Standard Pando uses the coefficient-table TF–target correlation.

Inspect the condition-GRN penalty fields before continuing:

```r
all_edges <- step1$grn_result$tf_peak_gene_condition_all

condition_active <- subset(
  all_edges,
  (is.na(analysis_mode) | analysis_mode == "condition_grn") &
    penalty_eligible %in% TRUE,
  select = c(
    tf, target, region, condition,
    estimate, padj, corr, corr_source,
    estimable, penalty_effect
  )
)

head(condition_active)
table(condition_active$corr_source, useNA = "ifany")
summary(abs(condition_active$corr))
summary(abs(condition_active$estimate))
```

When at least two cell-type jobs are available, Stage 1 distributes those jobs through `BPPARAM`. Every worker runs its own Pando job serially, so nested worker pools are not created. With one standard-Pando job, its existing target-level path may be used. A single condition-GRN job remains serial because pooled discovery, condition-specific discovery, dictionary freezing, and condition refits are treated as one coordinated contract.

```r
step1$grn_result$cell_type_analysis_mode
step1$grn_result$pando_execution_plan
step1$grn_result$pando_infer_argument_routing
```

## Stage 2: metacells

```r
step2 <- rc_regcompass_step_metacells(
  object = A,
  grn = step1,
  outdir = "run/02_metacells",
  condition_col = "condition",
  celltype_col = "cell_type",
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

Stage 2 reproduces the exact Stage 1 cell set. One WNN graph is built per broad cell type; final metacells remain condition-pure.

## Stage 3: biological reaction catalogue

```r
step3 <- rc_regcompass_step_meta_modules(
  grn = step1,
  metacells = step2,
  gem = gem,
  outdir = "run/03_meta_modules"
)
```

## Stage 4: Layer 1 reaction support

```r
step4 <- rc_regcompass_step_layer1(
  grn = step1,
  metacells = step2,
  meta_modules = step3,
  gem = gem,
  outdir = "run/04_layer1",
  gpr_and_method = "min",
  gene_half_saturation = 1
)
```

## Stage 5: Layer 2 metabolic scoring

```r
step5 <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "run/05_layer2",
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
  )
)
```

## Stage 6: result assembly

```r
result <- rc_regcompass_step_results(
  grn = step1,
  metacells = step2,
  meta_modules = step3,
  layer1 = step4,
  layer2 = step5,
  gem = gem,
  outdir = "run/06_results",
  species = "human"
)
```

To restart, load the last valid checkpoint and rerun only later stages. Do not combine checkpoints created from different cell sets, GEMs, media, or metadata columns.
