# Tutorial Level 2: stepwise run

Use this tutorial when each RegCompass stage should be run and saved independently.

## Stage 1: infer condition-by-cell-type GRNs

```r
step1 <- rc_regcompass_step_grn(
  object = A,
  gem = gem,
  outdir = "RegCompass_steps/01_grn",
  pfm = motif2tf,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  condition_col = "Group",
  celltype_col = "cell_type",
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
  parallel = TRUE,
  BPPARAM = upstream_bp
)

step1$grn_result$sample_status
```

## Stage 2: construct condition-level metacells

```r
step2 <- rc_regcompass_step_metacells(
  object = A,
  outdir = "RegCompass_steps/02_metacells",
  sample_col = NULL,
  condition_col = "Group",
  celltype_col = "cell_type",
  fragment_files = FALSE,
  metacell_args = list(
    gamma = 30,
    min_cells_per_stratum = 500,
    min_metacell_size = 10
  )
)

step2$pooled$metacell_meta
```

## Stage 3: construct reaction meta-modules

```r
step3 <- rc_regcompass_step_meta_modules(
  grn = step1,
  metacells = step2,
  gem = gem,
  outdir = "RegCompass_steps/03_meta_modules",
  layer1_args = list(
    top_k_neighbors = 5,
    min_shared_tfs = 1,
    min_tf_jaccard = 0,
    max_targets_per_tf = 200,
    expansion_mode = "ordered_once"
  )
)
```

Stage 3 projects metabolic-gene GRN components to complete-GPR core reactions, expands them through subsystem and direct KEGG/Reactome/master-Rhea annotations, and deduplicates reaction IDs across modules.

```r
catalogue <- step3$merged_modules
catalogue$merged_core_reactions
catalogue$merged_reaction_membership
table(catalogue$merged_reaction_membership$inclusion_stage)
```

## Stage 4: calculate integrated RNA+ATAC reaction support

```r
step4 <- rc_regcompass_step_layer1(
  metacells = step2,
  meta_modules = step3,
  gem = gem,
  outdir = "RegCompass_steps/04_layer1",
  regulatory_alpha = 1,
  tau = 0.20,
  parallel = TRUE,
  BPPARAM = upstream_bp
)
```

## Stage 5: build the medium-constrained model and score reactions

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
    time_limit = 600,
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

Stage 5 applies the selected medium, performs global FASTCORE completion, caches the model, and runs directional LP scoring. See [medium presets](medium-presets.md) for available presets and custom media.

```r
step5$model_cache_summary[, c(
  "medium_scenario",
  "n_merged_biological_reactions",
  "n_global_fastcore_support_reactions",
  "n_reactions",
  "file"
)]
```

## Stage 6: assemble annotated results

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

result$reaction_ranking
result$condition_contrast
result$merged_grn_meta_modules$merged_core_reactions
```
