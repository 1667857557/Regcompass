# Tutorial Level 2: true stepwise run with audit gates

Use this tutorial when each RegCompass stage must be run, inspected, saved, and restarted independently.

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

stopifnot(
  inherits(step1, "regcompass_grn_step"),
  all(step1$grn_result$sample_status$status == "ok")
)
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

stopifnot(
  inherits(step2, "regcompass_metacell_step"),
  setequal(
    colnames(step2$metacell_object),
    step2$pooled$metacell_meta$metacell_id
  )
)
```

## Stage 3: construct biological meta-modules

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

Stage 3 performs:

1. metabolic-gene GRN projection;
2. connected-component meta-module definition;
3. complete-GPR core mapping;
4. core-subsystem expansion;
5. KEGG/Reactome and master-Rhea reaction equivalence expansion;
6. reaction-ID deduplication across meta-modules.

It does **not** run FASTCORE and does **not** create a GEM.

```r
catalogue <- step3$merged_modules

stopifnot(
  identical(catalogue$is_gem, FALSE),
  identical(catalogue$fastcore_applied, FALSE),
  nrow(catalogue$merged_core_reactions) > 0,
  nrow(catalogue$merged_reaction_membership) > 0
)

table(catalogue$merged_reaction_membership$inclusion_stage)
```

The following old fields no longer exist:

```text
global_modules
global_core_reactions
global_reaction_membership
local_completed_reaction_membership
local_fastcore_summary
local_fastcore_diagnostics
```

## Stage 4: build integrated reaction support

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

stopifnot(inherits(step4, "regcompass_layer1_step"))
```

## Stage 5: build medium-specific union GEMs and score reactions

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

For each medium scenario, Stage 5:

1. applies the medium to the validated parent GEM;
2. starts from the merged biological reaction catalogue;
3. uses all merged complete-GPR reactions as targets;
4. performs the only FASTCORE completion;
5. saves one union GEM shared by all conditions and metacells;
6. runs directional two-step LP scoring.

```r
step5$model_cache_summary[, c(
  "medium_scenario",
  "n_merged_biological_reactions",
  "n_global_fastcore_support_reactions",
  "n_reactions",
  "build_strategy"
)]

union_gem <- readRDS(step5$model_cache_summary$file[[1]])
stopifnot(
  isTRUE(union_gem$is_union_gem),
  identical(
    union_gem$build_params$completion_stage,
    "single_global_fastcore_after_meta_module_merge"
  )
)
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

## Audit rule

A merged meta-module catalogue is not flux-completed and must not be called a union GEM. Only the cached, medium-constrained Stage 5 model is a union GEM.
