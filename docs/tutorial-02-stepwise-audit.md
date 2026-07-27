# Tutorial Level 2: stepwise run and audit

**Previous:** [Tutorial 1 — minimal one-shot run](tutorial-01-quick-start.md).

**Next:** [Tutorial 3 — restart and sensitivity](tutorial-03-advanced-restart.md).

This tutorial reproduces the one-shot workflow stage by stage and verifies that sample-aware bootstrap provenance reaches the expected Stage 1 outputs without changing later structural contracts.

## 1. Stage 1: shared-background GRN

```r
step1 <- rc_regcompass_step_grn(
  object = A,
  gem = gem,
  outdir = "RegCompass_steps/01_single_cell_grn",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  condition_col = "Group",
  celltype_col = "cell_type",
  sample_col = "sample_id",
  grn_mode = "multitask_shared_backbone",
  pando_args = list(
    min_cells = 100,
    pando_design_args = list(
      peak_to_gene_method = "GREAT",
      upstream = 100000,
      downstream = 0,
      extend = 1000000,
      only_tss = FALSE,
      min_tf_detection = 0,
      min_peak_detection = 0,
      min_target_detection = 0,
      max_edges_per_target = Inf
    )
  ),
  multitask_args = list(
    alpha = 0.5,
    nfolds = 5,
    n_bootstrap = 100,
    lambda_rule = "lambda.1se",
    min_selection_frequency = 0.7,
    min_sign_stability = 0.8,
    min_bootstrap_success_fraction = 0.8,
    candidate_screen_threshold = 0,
    max_edges_per_target = Inf,
    seed = 12345L
  ),
  parallel = TRUE,
  BPPARAM = upstream_bp,
  progress = TRUE
)
```

A valid `sample_col` activates condition-stratified sample-cluster bootstrap. Within each condition, donor/sample IDs are sampled with replacement, and all cells from every selected sample are retained. The model is re-centred within condition and refitted at the full-data selected lambda.

Audit the output contract:

```r
stopifnot(
  identical(step1$params$sample_col, "sample_id"),
  identical(
    step1$grn_result$bootstrap_policy$resampling_unit,
    "sample"
  )
)

with(step1$grn_result, {
  table(stability_diagnostics$bootstrap_resampling_unit)
  unique(stability_diagnostics$bootstrap_sample_col)
  celltype_fit_status[, c(
    "cell_type", "n_samples", "min_samples_per_condition",
    "bootstrap_resampling_unit", "bootstrap_fallback_reason"
  )]
  group_status[, c(
    "Group", "cell_type", "n_samples",
    "bootstrap_resampling_unit", "bootstrap_fallback_reason"
  )]
})
```

### Explicit fallback audit

```r
step1_fallback <- rc_regcompass_step_grn(
  object = A,
  gem = gem,
  outdir = "RegCompass_steps/01_single_cell_grn_fallback",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  condition_col = "Group",
  celltype_col = "cell_type",
  sample_col = "missing_sample_column",
  pando_args = list(min_cells = 100),
  multitask_args = list(n_bootstrap = 100),
  progress = TRUE
)
```

R prints a warning containing:

```text
Sample-aware bootstrap fallback: metadata column `missing_sample_column` does not exist. Falling back to condition-stratified cell resampling.
```

Then verify:

```r
stopifnot(
  identical(
    step1_fallback$grn_result$bootstrap_policy$resampling_unit,
    "cell"
  ),
  grepl(
    "does not exist",
    step1_fallback$grn_result$bootstrap_policy$fallback_reason,
    fixed = TRUE
  )
)
```

Do not use fallback for an existing column with missing/empty IDs. That input stops because silently dropping donors would invalidate cluster bootstrap.

## 2. Stage 2: condition-only metacells

```r
step2 <- rc_regcompass_step_metacells(
  object = A,
  outdir = "RegCompass_steps/02_condition_metacells",
  condition_col = "Group",
  celltype_col = "cell_type",
  fragment_files = FALSE,
  metacell_args = list(
    gamma = 30,
    rna_dims = 1:30,
    atac_dims = 2:30,
    min_cells_per_stratum = 500,
    min_metacell_size = 10,
    min_metacells_per_stratum = 2L,
    seed = 12345L,
    overwrite = FALSE
  ),
  progress = TRUE
)
```

`sample_col` is intentionally absent here. Stage 2 uses `strata_cols = condition_col`, passes `label = celltype_col`, and retains no artificial condition-pool metadata field.

## 3. Stage 3: condition cores and biological modules

```r
step3 <- rc_regcompass_step_meta_modules(
  grn = step1,
  metacells = step2,
  gem = gem,
  outdir = "RegCompass_steps/03_meta_modules",
  progress = TRUE
)
```

The sample-aware bootstrap changes which edges pass stability thresholds and can therefore change condition target genes and complete-GPR core reactions. It does not change the complete-GPR rule or annotation-expansion order.

## 4. Stage 4: integrated reaction support

```r
step4 <- rc_regcompass_step_layer1(
  metacells = step2,
  meta_modules = step3,
  gem = gem,
  outdir = "RegCompass_steps/04_layer1",
  regulatory_alpha = 1,
  gpr_and_method = "min",
  parallel = TRUE,
  BPPARAM = upstream_bp,
  progress = TRUE
)
```

## 5. Stage 5: shared union GEM and directional scoring

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
  BPPARAM = layer2_bp,
  progress = TRUE
)
```

All conditions and metacells still reuse the same shared medium-specific union GEM. Bootstrap mode affects regulatory evidence and core selection, not per-condition model reconstruction.

## 6. Stage 6: compact results

```r
result <- rc_regcompass_step_results(
  grn = step1,
  metacells = step2,
  meta_modules = step3,
  layer1 = step4,
  layer2 = step5,
  gem = gem,
  outdir = "RegCompass_steps/06_results",
  species = "human",
  progress = TRUE
)
```

`reaction_catalog`, `reaction_evidence`, `table_manifest`, and `stage_provenance` remain in the compact result; full coefficient and stability tables remain in detailed stage checkpoints.

Each stage prints its elapsed time in R after its final artifact is committed. Timing is not embedded in stage objects or written to timing TSV files. Metacell-level comparisons remain descriptive and are not biological-replicate treatment inference.
