# Tutorial 2: stepwise workflow and audit

Use the stepwise API to save, inspect, and restart each canonical stage. The
condition-GRN equations are in
[Tutorial 3](tutorial-03-mathematical-model.md); the public API is summarized in
[functions.md](functions.md).

## Progress and audit logs

All public stages accept `progress = TRUE` and persist `step_progress.tsv` and
`step_timing.tsv`. Stage 1 reports design resolution, normalization, motif
mapping, pooled and condition candidate discovery, exact edge union,
fixed-dictionary fitting, contract extraction, and artifact writing.

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
checkpointed target solver, nested CV, dense/matrix-free solver switch, or OOF
assignment. Errors stop immediately and the audit log records `stage_error`.

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

`BPPARAM = TRUE` is invalid. Do not place `parallel` or `BPPARAM` inside
`pando_infer_args`; the stage wrapper owns parallel execution. Layer 2 workers
force numerical libraries and HiGHS to one internal thread, so the worker count
controls process-level LP parallelism without nested thread multiplication.
Reduce `layer2_bp` workers first when memory is limiting.

## Stage 1: GRN inference

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

The Stage 1 threshold is fixed at 300 cells. A different supplied value is
overridden before normalization.

With at least two effective conditions, each retained cell type follows:

```text
pooled candidate discovery
+ each-condition candidate discovery
→ exact (target, TF, region) union
→ frozen target-specific dictionary
→ one unscaled Gaussian identity GLM per condition
→ complete coefficient and uncertainty table
→ estimable BH padj < 0.05 effects eligible for penalty
```

Candidate discovery uses Pando peak-to-gene domains, motif support,
peak-target correlation, and TF-target correlation. Candidate coefficients are
not reused as effects, and the final fit does not rerun correlation screening.

Inspect the fitted contracts and edge tables:

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
`estimable`, `zero_variance`, `aliased`, residual-df diagnostics, model `R²`, and
edge provenance. In multi-condition mode, `condition_active` is strictly the
estimable Pando edge set with BH `padj < 0.05`; no additional absolute-effect or
model-R² filter is applied to the penalty handoff. Those legacy extraction gates
remain relevant to the direct standard-Pando fallback.

Ordinary GLM P values are conditional on the frozen dictionary and do not
account for candidate-selection uncertainty.

With no condition column or fewer than two effective condition levels, Stage 1
uses `standard_pando` independently for each retained cell type and calculates
no condition coefficients.

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
  ),
  progress = TRUE,
  grn = step1
)
```

Passing `grn = step1` forces Stage 2 to reproduce the exact ordered fitted cell
set from Stage 1.

```r
step2$pooled$metacell_meta
step2$pooled$membership
step2$pooled$input_design
```

Each broad cell type has one independent multimodal WNN graph. Conditions share
that graph, and condition splits parent membership after clustering so final
metacells remain condition-pure.

## Stage 3: reaction catalogue and meta-modules

```r
step3 <- rc_regcompass_step_meta_modules(
  grn = step1,
  metacells = step2,
  gem = gem,
  outdir = "RegCompass_steps/03_meta_modules",
  progress = TRUE
)
```

```r
step3$condition_modules$core_gene_reaction
step3$merged_modules$merged_core_reactions
step3$merged_modules$merged_reaction_membership
```

Stage 3 derives condition-specific reaction evidence and stores the complete
condition module payload in a checksummed external RDS. The stage object retains
a lightweight reference so later stages do not duplicate the full payload in
memory.

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
  BPPARAM = upstream_bp,
  progress = TRUE
)
```

For multiple conditions, Stage 4 reconstructs `TF RNA × peak ATAC` on exact
paired cells, multiplies by the Pando `penalty_effect` from each cell's own
condition, sums by target, and then averages by exact SuperCell membership. It
does not refit, renormalize, or rebuild the dictionary at metacell level.

```r
step4$gene_regulatory_modifier
step4$reaction_expression
step4$projection_coverage
step4$projection_provenance
```

Historical output fields containing `_oof`, `common`, or `condition_unique` are
retained for API compatibility. The primary and common fields are aliases of the
BH-filtered fixed-dictionary projection; the condition-unique compatibility
matrix is zero.

## Stage 5: medium-specific directional penalties

Create one or more medium scenarios:

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
```

Custom reaction- or metabolite-level media described in
[Tutorial 1](tutorial-01-quick-start.md) can be passed in the same way.

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

All conditions and metacells use the same medium-specific structural model,
reaction order, bounds, target direction, target-flux fraction, and `vmax`.
Global FASTCORE completion is performed once per medium and cached.

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
  species = "human",
  progress = TRUE
)
```

```r
result$reaction_ranking
result$condition_contrast
result$microcompass$model_cache_summary
```

## Optional targeted reaction remapping

After Stage 5, selected reaction anchors can be mapped to direct KEGG,
Reactome, or master-Rhea equivalents and scored without rebuilding the cached
model. See
[Tutorial 4](tutorial-04-targeted-reaction-remapping.md).

## Restart rules

Changing candidate thresholds, target genes, conditions, cell types, RNA/ATAC
normalization, motif mapping, `padj_threshold`, rank handling, or the Stage 1
cell set invalidates Stage 1 and all downstream stages.

Changing metacell graph inputs or membership invalidates Stages 2–6 but does not
change already fitted single-cell GRNs.

Changing the reaction catalogue invalidates Stages 3–6.

Changing only medium scenarios or LP controls invalidates Stage 5 and final
results but does not require rerunning Stages 1–4.

Changing only the requested condition contrasts invalidates Stage 6 reporting,
not the GRN, metacell, or metabolic models.
