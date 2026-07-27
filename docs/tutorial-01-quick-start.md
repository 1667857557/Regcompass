# Tutorial Level 1: minimal one-shot run

This tutorial runs the complete RegCompassR 1.8.10 workflow with donor/sample-aware Stage 1 bootstrap.

**Next:** [Tutorial 2 — stepwise run and audit](tutorial-02-stepwise-audit.md).

## 1. Required metadata

The paired RNA+ATAC Seurat object requires complete condition and cell-type columns. A biological sample/donor column is optional but strongly recommended:

```text
condition_col   required
celltype_col    required
sample_col      optional; used only for Stage 1 bootstrap
```

When `sample_col` names a valid metadata column, RegCompass resamples sample/donor clusters with replacement separately inside each condition. Every selected sample contributes all of its cells. The number of cells in a bootstrap replicate may therefore vary with sample size.

When `sample_col = NULL`, or the named column does not exist, Stage 1 prints the exact reason and falls back to condition-stratified cell resampling. An existing sample column containing `NA` or empty identifiers is an error. Stage 2 remains condition-only and does not use `sample_col` for metacell construction.

## 2. Prepare the GEM and medium

```r
library(RegCompassR)
library(Seurat)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)

gem <- rc_prepare_gem(
  species = "human",
  version = "2.0.0",
  source = "bundled"
)

medium_scenarios <- rc_make_medium_scenarios(
  gem,
  scenario = "high_glucose",
  species = "human"
)
```

## 3. Run the complete workflow

```r
result <- rc_run_regcompass_one_shot(
  object = A,
  outdir = "RegCompass_result",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  gem = gem,
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
    global_penalty_factor = 1,
    deviation_penalty_factor = 1,
    lambda_rule = "lambda.1se",
    nfolds = 5,
    n_bootstrap = 100,
    min_selection_frequency = 0.7,
    min_sign_stability = 0.8,
    min_bootstrap_success_fraction = 0.8,
    min_cv_rsq = 0,
    candidate_screen_threshold = 0,
    max_edges_per_target = Inf,
    min_detected_cells_per_condition = 10,
    min_detection_fraction_per_condition = 0.01,
    seed = 12345L
  ),

  fragment_files = FALSE,
  metacell_args = list(
    rna_reduction = "pca",
    rna_dims = 1:30,
    atac_reduction = "lsi",
    atac_dims = 2:30,
    gamma = 30,
    seed = 12345L,
    min_cells_per_stratum = 500,
    min_metacell_size = 10,
    min_metacells_per_stratum = 2L,
    overwrite = FALSE
  ),

  layer1_args = list(
    regulatory_alpha = 1,
    gpr_and_method = "min"
  ),
  medium_scenarios = medium_scenarios,
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
  upstream_workers = 6,
  layer2_workers = 30,
  progress = TRUE
)
```

Each stage prints a line such as:

```text
RegCompass timing: single_cell_grn [success] 00:12:34.567
```

Timing is no longer stored in `result$timing`, `step_timing.tsv`, or `00_execution_timing.tsv`.

## 4. Audit the bootstrap policy

```r
step1 <- readRDS(
  "RegCompass_result/01_single_cell_grn/step_grn.rds"
)

step1$params$sample_col
step1$grn_result$bootstrap_policy
head(step1$grn_result$stability_diagnostics[, c(
  "bootstrap_method",
  "bootstrap_resampling_unit",
  "bootstrap_sample_col",
  "n_bootstrap_samples_total",
  "min_bootstrap_samples_per_condition",
  "bootstrap_fallback_reason"
)])
```

Expected sample-aware method:

```text
condition_stratified_sample_cluster_nonparametric
```

Fallback method:

```text
condition_stratified_cell_nonparametric_fallback
```

The same provenance fields are also available in `celltype_fit_status`, `group_status`, and `bootstrap_stability_diagnostics.tsv.gz`.

## 5. Interpret the workflow boundary

`sample_col` changes only the bootstrap estimate of edge reproducibility. It does not change:

- the shared Pando candidate universe;
- condition-centred direct-theta fitting;
- Stage 2 condition-only metacells;
- complete-GPR condition core definitions;
- the merged reaction catalogue;
- the single shared medium-specific union GEM reused across all conditions and metacells.

The final object remains compact through `table_manifest` and `stage_provenance`; full bootstrap and model diagnostics remain in detailed stage checkpoints. Metacell-level comparisons are not biological-replicate treatment inference.

Continue to [Tutorial 2](tutorial-02-stepwise-audit.md) for the stage-by-stage workflow.
