# Tutorial 2: stepwise workflow

Use the stepwise API to save, inspect and restart each canonical stage. The
condition-GRN contract is in
[condition-comparable-grn.md](condition-comparable-grn.md). Public API:
[functions.md](functions.md).

## Progress and audit log

All public stages accept `progress = TRUE` and persist `step_progress.tsv` and
`step_timing.tsv`. Stage 1 reports design resolution, global normalization,
candidate initialization, motif mapping, global/condition candidate discovery,
exact edge union, fixed-dictionary GLM fitting, contract extraction and artifact
writing.

```r
options(RegCompassR.progress = TRUE)
step1 <- rc_regcompass_step_grn(..., progress = TRUE)
progress_log <- read.delim(
  "RegCompass_steps/01_grn/step_progress.tsv",
  check.names = FALSE
)
progress_log[, c("phase", "elapsed_hms", "detail", "context")]
```

The current condition route has no native condition runtime, lambda path,
checkpointing, nested CV, dense/matrix-free solver choice or OOF assignment.
Errors stop immediately and the final audit row records `stage_error`.

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

`BPPARAM = TRUE` is invalid. Do not put `parallel` or `BPPARAM` inside
`pando_infer_args`.

## Stage 1: common-dictionary condition GRNs

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
    min_cells = 300L,
    min_abs_estimate = 0,
    min_model_rsq = 0.1,
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

The fixed Stage 1 threshold is 300 cells. A different supplied value is
overridden by the hardening contract before normalization.

For each retained cell type, Stage 1 runs:

```text
global candidate discovery
+ each-condition candidate discovery
→ exact (target, TF, region) union
→ frozen target-specific dictionary
→ one unscaled Gaussian identity GLM per condition
→ BH-adjusted coefficient table
→ padj < 0.05 effects eligible for penalty
```

Candidate discovery uses Pando peak-to-gene domains, motif support,
peak-target correlation and TF-target correlation. The first-stage candidate
coefficients are not reused. The final fit never re-runs correlation screening.

Inspect the complete and penalty-eligible effects:

```r
step1$params$analysis_mode
step1$grn_result$condition_grn_fits
step1$grn_result$tf_peak_gene_condition_effect_all
step1$grn_result$tf_peak_gene_condition_effect
step1$grn_result$condition_fit_status
```

Stage 1 artifacts include:

```text
pando_group_status.tsv.gz
pando_tf_peak_gene_condition_all.tsv.gz
pando_tf_peak_gene_condition_active.tsv.gz
pando_tf_peak_gene_universal.tsv.gz
pando_condition_grn_fits.rds
```

`condition_all` retains `estimate`, `std_err`, `statistic`, `pval`, `padj`,
`estimable`, `zero_variance`, `aliased`, model `R²`, and the exact edge
provenance. `condition_active` requires BH significance, estimability, the
configured absolute-effect gate and `min_model_rsq`.

Ordinary GLM P values are conditional on the frozen dictionary; they do not
account for candidate-selection uncertainty.

With no condition column or fewer than two condition levels, Stage 1 uses
`standard_pando` independently for each retained cell type and calculates no
condition coefficients.

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
    k.knn = 30L,
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

Each broad cell type has one independent multimodal WNN graph. Conditions share
that graph, and condition splits the parent membership after clustering so final
metacells remain condition-pure.

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

## Stage 4: regulatory reaction support

```r
step4 <- rc_regcompass_step_layer1(
  grn = step1,
  metacells = step2,
  meta_modules = step3,
  gem = gem,
  outdir = "RegCompass_steps/04_layer1",
  projection_component = "condition",
  regulatory_alpha = 1,
  gpr_and_method = "min",
  parallel = TRUE,
  BPPARAM = upstream_bp
)
```

For multiple conditions, Stage 4 reconstructs `TF RNA × peak ATAC` on the exact
paired cells, multiplies by Pando `penalty_effect`, sums by target, and then
aggregates using the exact SuperCell membership table. No coefficient,
normalization or edge dictionary is recomputed at metacell level.

```r
step4$gene_regulatory_modifier
step4$reaction_expression
step4$projection_coverage
step4$projection_provenance
```

Historical output fields containing `_oof`, `common`, or `condition_unique` are
retained for API compatibility. In the current method, the primary and common
fields are aliases of the BH-filtered fixed-dictionary full-fit projection, and
the condition-unique compatibility decomposition is zero.

## Stage 5: media and directional penalty

```r
medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = c(
    "normal_human_plasma",
    "high_glucose",
    "low_glucose"
  ),
  species = "human"
)

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

All conditions and metacells use the same medium-specific model, reaction order,
bounds, target directions and `vmax`.

```r
step5$penalty
step5$vmax
step5$model_cache_summary
```

Historical `penalty_condition_full_oof` and `penalty_common_oof` fields are
compatibility aliases of the current primary penalty; the condition-unique
increment is zero.

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
```

## Restart rules

Changing candidate thresholds, target genes, conditions, cell types, RNA/ATAC
normalization, motif mapping, `padj_threshold`, rank handling or Stage 1 cell
sets invalidates Stage 1 and every downstream stage.

Changing metacell graph inputs or membership invalidates Stages 2-6 but does not
change the already fitted single-cell GRNs.

Changing only medium scenarios or LP controls invalidates Stage 5 and final
results but does not require rerunning Stages 1-4.
