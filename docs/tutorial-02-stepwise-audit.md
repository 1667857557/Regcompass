# Tutorial Level 2: stepwise run and scientific audit

**Previous:** [Tutorial 1 — minimal one-shot run](tutorial-01-quick-start.md).

**This tutorial:** reproduces the same RegCompassR 1.8.9 workflow as six explicit stages. Use it when intermediate tables, mathematical invariants, or restart points must be inspected.

**Next:** [Tutorial 3 — restart and sensitivity](tutorial-03-advanced-restart.md) reuses these stage objects. [Tutorial 5](tutorial-05-condition-differential-analysis.md) connects the final tables back to Stage 1–5 evidence.

Use stable object names throughout:

```text
step1 = GRN
step2 = metacells
step3 = meta-modules
step4 = Layer 1 evidence
step5 = Layer 2 shared-model scoring
result = compact Stage 6 output
```

## 1. Parallel backends

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

Stage 1 parallelises cell types; target-level `glmnet` fits remain single-threaded. Stage 4 can reuse `upstream_bp`; Stage 5 uses `layer2_bp`.

## 2. Stage 1 — shared structural GRN and condition sub-GRNs

```r
step1 <- rc_regcompass_step_grn(
  object = A,
  gem = gem,
  outdir = "RegCompass_steps/01_grn",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  condition_col = "Group",
  celltype_col = "cell_type",
  grn_mode = "multitask_shared_backbone",
  pando_args = list(
    min_cells = 100,
    pando_design_args = list(
      peak_to_gene_method = "Signac",
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
  parallel = TRUE,
  BPPARAM = upstream_bp
)
```

Detailed Stage 1 tables:

```r
step1$grn_result$celltype_fit_status
step1$grn_result$group_status
step1$grn_result$tf_peak_gene_candidates
step1$grn_result$tf_peak_gene_global
step1$grn_result$tf_peak_gene_condition_all
step1$grn_result$tf_peak_gene_significant
step1$grn_result$condition_target_genes
step1$grn_result$target_model_diagnostics
step1$grn_result$stability_diagnostics
```

The candidate table distinguishes two universes:

```text
edge_universe_id       = Pando structural motif/domain universe
model_edge_universe_id = condition-aware observable model universe
model_observable       = structural edge retained for fitting
```

For condition \(c\), the observability threshold is

\[
m_c=\min\left(n_c,\max\left(10,\left\lceil0.01n_c\right\rceil\right)\right).
\]

An edge enters the shared model universe when the non-zero TF-RNA × peak-ATAC predictor and target RNA each meet \(m_c\) in at least one condition. This uses detection only; it does not use target correlation or effect size.

For each cell type, every condition must use the same `model_edge_universe_id`. The fitted fields satisfy:

```text
effective_estimate = global_estimate + condition_deviation
stable_estimate = effective_estimate × selection_frequency × sign_stability
```

and for every edge:

\[
\sum_c\delta_{e,c}=0.
\]

Global and condition-deviation coordinates use equal explicit penalty factors by default. Values above one are explicit sensitivity priors for a more conserved shared backbone.

The bootstrap is full-size and condition stratified:

```text
for each bootstrap and each condition:
  sample n_c cells with replacement
  recalculate condition centres
  apply the shared edge scale
  fit at the full-data selected lambda
```

Check predictive validity and bootstrap completion before interpreting stability:

```r
with(
  step1$grn_result$target_model_diagnostics,
  summary(cv_rsq)
)

with(
  step1$grn_result$target_model_diagnostics,
  summary(n_bootstrap_success / n_bootstrap_requested)
)

step1$grn_result$stability_diagnostics[, intersect(c(
  "selection_frequency",
  "selection_frequency_mc_se",
  "selection_frequency_lower_95",
  "selection_frequency_upper_95",
  "sign_stability",
  "sign_agreement_fraction",
  "cv_predictive_above_null",
  "active_edge"
), colnames(step1$grn_result$stability_diagnostics)), drop = FALSE]
```

An active edge must have strictly positive out-of-fold \(R^2\). `min_sign_stability = 0.8` corresponds to at least 90% majority-sign agreement among selected bootstrap fits because \(\rho=|2q-1|\).

See [Pando and multitask GRN parameter policy](grn-parameter-policy.md) for the parameter derivations and supported sensitivity analyses.

Stage 1 checkpoint:

```text
RegCompass_steps/01_grn/step_grn.rds
```

## 3. Stage 2 — condition strata plus exact SuperCell2 cell-type labels

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
    gamma = 30,
    seed = 12345L,
    min_cells_per_stratum = 500,
    min_metacell_size = 10,
    min_metacells_per_stratum = 2L,
    overwrite = FALSE
  )
)
```

The current contract is:

```text
RegCompass strata_cols = condition_col
SuperCell2 label = celltype_col
```

There is no `sample_col`, artificial condition-pool column, sample balancing, or sample-level grouping.

```r
formals(rc_make_supercell2_metacells)
step2$pooled$strata_cols
step2$pooled$label_col
step2$pooled$metacell_meta
step2$pooled$celltype_composition_summary
```

Every metacell must be label pure:

```r
stopifnot(
  all(step2$pooled$celltype_composition_summary$n_celltypes == 1L),
  !any(step2$pooled$celltype_composition_summary$mixed_celltype_metacell)
)
```

Stage 2 checkpoint:

```text
RegCompass_steps/02_metacells/step_metacells.rds
RegCompass_steps/02_metacells/merged_metacell_object.rds
```

## 4. Stage 3 — condition-specific complete-GPR cores

```r
step3 <- rc_regcompass_step_meta_modules(
  grn = step1,
  metacells = step2,
  gem = gem,
  outdir = "RegCompass_steps/03_meta_modules"
)
```

For each `condition × cell type` `group_id`:

```text
bootstrap-active edges
→ regulated metabolic targets
→ complete GPR branches
→ core reactions
→ core subsystems
→ direct KEGG/Reactome equivalents
→ direct master-Rhea equivalents
→ biological meta-module
```

A reaction is a core when at least one branch is complete. The branch table distinguishes:

```text
reaction_is_core = any complete branch exists
is_core/group_complete = this branch is complete
```

```r
step3$condition_modules$supported_metabolic_genes
step3$condition_modules$core_gene_reaction
step3$condition_modules$reaction_membership
step3$merged_modules$merged_core_reactions
step3$merged_modules$merged_reaction_membership
```

The merged output is a biological reaction catalogue, not yet a GEM.

Stage 3 checkpoint:

```text
RegCompass_steps/03_meta_modules/step_meta_modules.rds
```

## 5. Stage 4 — RNA support, ATAC modifier, and GPR aggregation

```r
step4 <- rc_regcompass_step_layer1(
  metacells = step2,
  meta_modules = step3,
  gem = gem,
  outdir = "RegCompass_steps/04_layer1",
  regulatory_alpha = 1,
  gpr_and_method = "min",
  parallel = TRUE,
  BPPARAM = upstream_bp
)
```

```r
step4$gene_support_rna
step4$gene_regulatory_modifier
step4$gene_support_multiome
step4$reaction_expression
step4$gpr_diagnostics
```

TFs sharing one measured peak are signed-summed before ATAC projection. The denominator is shared across conditions for each target. When no active edge exists, the modifier is zero and multiome support equals RNA support exactly.

The canonical GPR rule is:

```text
AND = min
OR = additive isozyme support
```

Stage 4 checkpoint:

```text
RegCompass_steps/04_layer1/step_layer1.rds
```

## 6. Stage 5 — one shared union GEM per medium

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

Stage 5 performs the only FASTCORE completion. Within each medium, all conditions and metacells reuse identical reaction IDs, `S`, `lb`, and `ub`.

```r
step5$model_cache_summary
step5$source_core_reactions
step5$source_merged_reaction_membership
step5$union_gem_policy
```

`completion_time_limit` applies only to union-GEM construction, not to directional scoring LPs.

Stage 5 checkpoint:

```text
RegCompass_steps/05_layer2/step_layer2.rds
```

## 7. Stage 6 — assemble a compact final result

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
result$table_manifest
result$reaction_ranking
result$condition_contrast
result$active_regulatory_edges
result$condition_target_genes
result$core_reactions
result$reaction_catalog
result$reaction_evidence
```

The final object does not duplicate the full Stage 1–4 objects:

```r
result$stage_provenance$detailed_intermediates_embedded
result$stage_provenance$detailed_sources
```

Use `result` for ranking, condition analysis, selection, and plotting. Use `step1`–`step5` when auditing all candidates, all coefficients, matrices, or complete membership tables.

## 8. Handoff to later tutorials

For sensitivity or restart, continue to [Tutorial 3](tutorial-03-advanced-restart.md) with `step1`–`step5`.

For direct database-equivalent non-core targets, continue to [Tutorial 4](tutorial-04-targeted-reaction-remapping.md) with `step3`, `step4`, and `step5`.

For biological interpretation and condition comparison, continue to [Tutorial 5](tutorial-05-condition-differential-analysis.md) with `result` and optionally the stage objects.
