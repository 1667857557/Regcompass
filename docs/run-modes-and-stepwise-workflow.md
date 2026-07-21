# RegCompassR execution modes

RegCompassR supports two equivalent execution modes:

1. **One-shot mode** runs the canonical workflow through one function call.
2. **Stepwise mode** exposes five public stages so that each output can be inspected, saved, adjusted and selectively rerun.

Both modes use the same condition-pooled biological design, Pando/meta-module logic, Layer 1 evidence model, shared structural GEM and directional COMPASS-like LP. Stepwise mode changes execution control, not the mathematical model.

## Shared setup

```r
library(RegCompassR)
library(Pando)
library(BSgenome.Hsapiens.UCSC.hg38)

data(motifs, package = "Pando")
data(SCREEN.ccRE.UCSC.hg38, package = "Pando")

gem <- rc_prepare_gem(
  species = "human",
  version = "2.0.0"
)

medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "normal_human_plasma",
  species = "human"
)

bp <- BiocParallel::SnowParam(workers = 8, type = "SOCK")
```

The input Seurat object must contain RNA and ATAC assays and the metadata columns `sample_id`, `condition` and `cell_type`. Each biological sample must map to exactly one condition. The package warns when fewer than two biological samples are available in a condition because metacells are descriptive pseudo-observations rather than biological replicates.

## Mode A: one-shot execution

Use one-shot mode for routine analyses after the parameters have been established.

```r
result <- rc_run_regcompass_one_shot(
  object = object,
  outdir = "RegCompass_result",
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  fragment_files = FALSE,
  species = "human",
  gem = gem,
  medium_scenarios = medium_scenarios,
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type",
  model_mode = "meta_module_gem",
  metacell_args = list(
    gamma = 20,
    min_cells_per_stratum = 500,
    min_metacell_size = 10
  ),
  pando_args = list(
    min_metacells = 10,
    pando_initiate_args = list(
      regions = SCREEN.ccRE.UCSC.hg38
    ),
    pando_infer_args = list(
      method = "glm",
      tf_cor = 0.1,
      peak_cor = 0,
      adjust_method = "fdr"
    )
  ),
  layer1_args = list(
    regulatory_alpha = 1,
    tau = 0.20,
    local_fastcore = TRUE
  ),
  layer2_args = list(
    target_direction = "both",
    solver = "highs"
  ),
  upstream_workers = 8,
  layer2_workers = 8
)
```

The final object is also written to `RegCompass_result/regcompass_result.rds`.

## Mode B: stepwise execution

Use stepwise mode while choosing parameters, diagnosing data problems or auditing intermediate biological results.

### Step 1: condition-pooled metacells

```r
step1 <- rc_regcompass_step_metacells(
  object = object,
  outdir = "RegCompass_steps/01_metacells",
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type",
  rna_assay = "RNA",
  atac_assay = "ATAC",
  fragment_files = FALSE,
  metacell_args = list(
    gamma = 20,
    min_cells_per_stratum = 500,
    min_metacell_size = 10,
    BPPARAM = bp
  )
)
```

Inspect before continuing:

```r
dim(step1$metacell_object)
head(step1$pooled$metacell_meta)
head(step1$pooled$sample_composition)
with(step1$pooled$metacell_meta, table(condition, cell_type))
summary(step1$pooled$metacell_meta$effective_sample_n)
step1$metacell_object@misc$regcompass_atac_normalization
stopifnot(
  setequal(
    colnames(step1$metacell_object),
    step1$pooled$metacell_meta$metacell_id
  )
)
```

Peaks with zero total counts across the pooled metacell object are removed before the cell-type-shared TF-IDF calculation. Inspect `n_zero_count_peaks_excluded` and `n_retained_peaks` in the normalization metadata.

Adjust `gamma`, `min_cells_per_stratum` or `min_metacell_size` here. Any change to Step 1 invalidates Steps 2–5.

Checkpoint: `RegCompass_steps/01_metacells/step_metacells.rds`.

### Step 2: Pando GRNs and meta-modules

```r
step2 <- rc_regcompass_step_meta_modules(
  metacells = step1,
  gem = gem,
  outdir = "RegCompass_steps/02_meta_modules",
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pando_args = list(
    min_metacells = 10,
    pando_initiate_args = list(
      regions = SCREEN.ccRE.UCSC.hg38
    ),
    pando_infer_args = list(
      method = "glm",
      tf_cor = 0.1,
      peak_cor = 0,
      adjust_method = "fdr"
    ),
    padj_threshold = 0.05,
    min_model_rsq = 0.1
  ),
  layer1_args = list(
    local_fastcore = TRUE,
    local_fastcore_args = list(
      solver = "highs"
    )
  ),
  parallel = TRUE,
  BPPARAM = bp
)
```

Inspect before continuing:

```r
with(step2$condition_modules$sample_status, table(status))
step2$condition_modules$sample_status[, c(
  "group_id", "n_atac_peaks_input", "n_zero_count_peaks_excluded",
  "n_atac_peaks_used"
)]
head(step2$condition_modules$tf_peak_gene_significant)
summary(step2$condition_modules$tf_peak_gene_significant$rsq)
head(step2$condition_modules$core_gene_reaction)
head(step2$condition_modules$local_fastcore_summary)
length(unique(step2$global_modules$global_core_reactions$reaction_id))
length(unique(step2$global_modules$global_reaction_membership$reaction_id))
```

Every condition-by-cell-type Pando fit must have status `ok`. Peaks with zero counts within a Pando group are removed before motif and GRN inference even when they were nonzero in another condition. Check that retained edges have finite `rsq`, that complete-GPR core reactions exist, and that local FASTCORE support is not unexpectedly large.

Adjust Pando thresholds, motif/region inputs or local FASTCORE settings here. Any change to Step 2 invalidates Steps 3–5 but does not require rebuilding metacells.

Checkpoint: `RegCompass_steps/02_meta_modules/step_meta_modules.rds`.

### Step 3: integrated Layer 1 reaction expression

```r
step3 <- rc_regcompass_step_layer1(
  metacells = step1,
  meta_modules = step2,
  gem = gem,
  outdir = "RegCompass_steps/03_layer1",
  regulatory_alpha = 1,
  tau = 0.20,
  gene_half_saturation = 1,
  parallel = TRUE,
  BPPARAM = bp
)
```

Inspect before continuing:

```r
dim(step3$gene_support_rna)
dim(step3$gene_regulatory_modifier)
dim(step3$gene_support_multiome)
dim(step3$reaction_expression)
range(step3$gene_support_rna, na.rm = TRUE)
range(step3$gene_regulatory_modifier, na.rm = TRUE)
range(step3$gene_support_multiome, na.rm = TRUE)
summary(as.numeric(step3$reaction_expression))
with(step3$gpr_diagnostics, table(capacity_missing_flag))
stopifnot(
  identical(
    colnames(step3$reaction_expression),
    as.character(step3$unit_meta$pool_id)
  )
)
```

Expected ranges are `[0,1]` for RNA and multiome gene support and `[-1,1]` for the regulatory modifier. Reaction expression is non-negative and is not a flux. Layer 1 parallelizes reaction-level GPR aggregation.

Adjust `regulatory_alpha`, `tau` or `gene_half_saturation` here. Any change to Step 3 invalidates Steps 4–5 only.

Checkpoint: `RegCompass_steps/03_layer1/step_layer1.rds`.

### Step 4: shared-GEM directional LP scoring

```r
step4 <- rc_regcompass_step_layer2(
  layer1 = step3,
  meta_modules = step2,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "RegCompass_steps/04_layer2",
  model_mode = "meta_module_gem",
  layer2_args = list(
    target_direction = "both",
    omega = 0.95,
    solver = "highs",
    time_limit = 60
  ),
  parallel = TRUE,
  BPPARAM = bp
)
```

Inspect before finalizing:

```r
dim(step4$penalty)
dim(step4$vmax)
mean(step4$feasible)
with(step4$lp_diagnostics, table(step2_status, useNA = "ifany"))
head(step4$model_cache_summary)
head(step4$model_diagnostics)
with(step4$penalty_components, table(missing_expression_flag))
stopifnot(
  identical(colnames(step4$penalty), colnames(step3$reaction_expression))
)
```

For expression-linked reactions, non-finite reaction expression is zero-filled before penalty calculation. Unmeasured expression and explicit zero expression therefore receive the same strictest expression-linked cost of 1. Missingness remains separately recorded in `penalty_components$missing_expression_flag`.

Review blocked targets, solver failures, medium mapping and the number of biological versus FASTCORE support reactions. Adjust `model_mode`, medium constraints, `omega`, target direction, solver or time limit here. Any change to Step 4 invalidates Step 5 only.

Checkpoint: `RegCompass_steps/04_layer2/step_layer2.rds`.

### Step 5: ranking and canonical result assembly

```r
result <- rc_regcompass_step_results(
  metacells = step1,
  meta_modules = step2,
  layer1 = step3,
  layer2 = step4,
  gem = gem,
  outdir = "RegCompass_steps/final",
  species = "human"
)
```

Inspect final outputs:

```r
head(result$reaction_ranking)
head(result$condition_summary)
head(result$condition_contrast)
table(result$reaction_ranking$cell_type, result$reaction_ranking$condition)
```

Checkpoint: `RegCompass_steps/final/regcompass_result.rds`.

## Parallel controls

Steps 2–4 accept both `parallel` and `BPPARAM`:

- `parallel = FALSE` forces sequential execution;
- `parallel = TRUE, BPPARAM = NULL` uses the default RegCompass backend;
- `BPPARAM = FALSE` forces base `lapply()`;
- a `BiocParallelParam` object selects an explicit backend.

Logical `TRUE` is not a valid `BPPARAM` object. Use `BiocParallel::SnowParam()`, `BiocParallel::MulticoreParam()` or `FALSE`.

## Restarting from checkpoints

Each stage writes a complete RDS object. A later session can restart from the latest valid stage:

```r
step1 <- readRDS("RegCompass_steps/01_metacells/step_metacells.rds")
step2 <- readRDS("RegCompass_steps/02_meta_modules/step_meta_modules.rds")
step3 <- readRDS("RegCompass_steps/03_layer1/step_layer1.rds")
```

Then rerun Step 4 with a different medium, solver or structural model without repeating Pando inference.

## Dependency rule for parameter changes

| Changed parameter or input | First stage to rerun | Downstream stages to rerun |
|---|---:|---|
| Cells, metadata, assays, metacell settings | 1 | 2–5 |
| PFM, genome, Pando settings, GRN thresholds, local FASTCORE | 2 | 3–5 |
| `regulatory_alpha`, `tau`, RNA half-saturation | 3 | 4–5 |
| Medium, structural mode, target direction, `omega`, solver | 4 | 5 |
| Ranking display or downstream plotting only | 5 or external analysis | none |

Do not combine checkpoints produced from different upstream parameterizations. The stage functions preserve metadata-column configuration and reject incompatible metacell and meta-module stage objects.
