# Tutorial Level 3: restart, sensitivity, and diagnostics

**Previous:** [Tutorial 2 — stepwise run and audit](tutorial-02-stepwise-audit.md) creates `step1`–`step5` and the compact `result` used here.

**This tutorial:** identifies the earliest stage that must be rerun after a parameter or data change, then rebuilds every dependent stage.

**Next:** [Tutorial 4](tutorial-04-targeted-reaction-remapping.md) reuses an unchanged Stage 5 union GEM for selected non-core reactions. [Tutorial 5](tutorial-05-condition-differential-analysis.md) interprets the rebuilt compact result.

## 1. Load a completed stepwise run

```r
step1 <- readRDS("RegCompass_steps/01_grn/step_grn.rds")
step2 <- readRDS("RegCompass_steps/02_metacells/step_metacells.rds")
step3 <- readRDS("RegCompass_steps/03_meta_modules/step_meta_modules.rds")
step4 <- readRDS("RegCompass_steps/04_layer1/step_layer1.rds")
step5 <- readRDS("RegCompass_steps/05_layer2/step_layer2.rds")
result <- readRDS("RegCompass_steps/06_results/regcompass_result.rds")
```

The compact result records where detailed information is stored:

```r
result$stage_provenance$detailed_sources
```

Never combine stages from different cells, metadata, GEM fingerprints, or assay states. Stage validators intentionally reject incompatible objects.

## 2. Restart dependency graph

```text
Stage 1 GRN ───────┐
                    ├→ Stage 3 modules → Stage 4 evidence → Stage 5 scoring → Stage 6 result
Stage 2 metacells ─┘
```

### Rerun Stage 1 and Stage 3–6 after changing

- motifs, genome, regulatory regions, or peak-to-gene domains;
- Pando structural detection thresholds;
- condition-aware observability thresholds;
- direct-theta elastic-net settings, CV rule, bootstrap number, or active-edge thresholds;
- condition/cell-type metadata used by Stage 1;
- single-cell RNA or ATAC matrices;
- GEM GPR target genes.

Stage 2 can be reused only when cells, labels, assays, reductions, and metacell parameters are unchanged.

### Rerun Stage 2 and Stage 3–6 after changing

- cells or RNA/ATAC count matrices;
- condition or cell-type labels;
- PCA/LSI reductions or dimensions;
- `gamma`, seed, minimum cell count, or minimum metacell size;
- fragment inputs.

The current Stage 2 contract has no sample column:

```text
RegCompass splits by condition
SuperCell2 receives label = celltype_col
```

### Rerun Stage 3–6 after changing

- Stage 1 active edges or condition target genes;
- GPR rules;
- subsystem mapping;
- KEGG, Reactome, or master-Rhea cross-references.

### Rerun Stage 4–6 after changing

- `regulatory_alpha`;
- `gpr_and_method`;
- RNA/ATAC half-saturation settings;
- metacell RNA or ATAC evidence;
- stable GRN projection weights.

### Rerun Stage 5–6 after changing

- medium composition or exchange bounds;
- target direction, `omega`, solver, or flux threshold;
- union-GEM completion parameters;
- Stage 3 reaction membership;
- Stage 4 reaction penalties.

### Rerun only Stage 6 after changing

- compact output schema or exported table policy;
- reaction annotation presentation;
- final output directory.

Stage 6 does not refit biology or LP models.

## 3. GRN sensitivity with an unchanged GREAT structural universe

The following analysis changes regularisation, observability, bootstrap size and
activation thresholds while preserving the canonical GREAT structural candidate
definition.

```r
step1_sensitive <- rc_regcompass_step_grn(
  object = A,
  gem = gem,
  outdir = "RegCompass_restart/01_grn_sensitive",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  condition_col = "Group",
  celltype_col = "cell_type",
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
    alpha = 0.75,
    global_penalty_factor = 1,
    deviation_penalty_factor = 1,
    lambda_rule = "lambda.1se",
    nfolds = 5,
    n_bootstrap = 200,
    min_selection_frequency = 0.8,
    min_sign_stability = 0.9,
    min_bootstrap_success_fraction = 0.8,
    min_cv_rsq = 0.05,
    candidate_screen_threshold = 0,
    max_edges_per_target = Inf,
    min_detected_cells_per_condition = 20,
    min_detection_fraction_per_condition = 0.02,
    seed = 12345L
  )
)
```

This is deliberately more selective than the canonical model:

- `alpha = 0.75` increases the L1 contribution acting directly on condition-specific \(\theta_{e,c}\);
- both compatibility penalty aliases remain equal because there is one direct-theta penalty;
- `min_cv_rsq = 0.05` requires stronger out-of-fold predictive value;
- the observability rule changes from `max(10, 1%)` to `max(20, 2%)`;
- 200 bootstrap fits reduce Monte Carlo error but do not add biological replicates.

Unequal settings such as `global_penalty_factor = 1` and
`deviation_penalty_factor = 3` are rejected. The direct-theta model has no
separate latent deviation-coordinate penalty; a conserved-backbone prior would
require a separate fused or grouped estimator.

Compare both universe identifiers:

```r
structural_old <- unique(
  step1$grn_result$tf_peak_gene_candidates$edge_universe_id
)
structural_new <- unique(
  step1_sensitive$grn_result$tf_peak_gene_candidates$edge_universe_id
)

model_old <- unique(na.omit(
  step1$grn_result$tf_peak_gene_candidates$model_edge_universe_id
))
model_new <- unique(na.omit(
  step1_sensitive$grn_result$tf_peak_gene_candidates$model_edge_universe_id
))

structural_old
structural_new
model_old
model_new
```

Expected interpretation:

```text
structural_old == structural_new
model_old may differ from model_new
```

The Pando structural universe should remain unchanged because its GREAT domain,
motif, and pooled detection settings are unchanged. The model universe may
change because the condition-aware observability threshold changed.

Check predictive validity and bootstrap completion:

```r
step1_sensitive$grn_result$target_model_diagnostics[, intersect(c(
  "target",
  "cv_rsq",
  "cv_predictive_above_null",
  "n_bootstrap_requested",
  "n_bootstrap_success",
  "bootstrap_success_fraction",
  "n_active_condition_edges"
), colnames(step1_sensitive$grn_result$target_model_diagnostics)), drop = FALSE]

step1_sensitive$grn_result$stability_diagnostics[, intersect(c(
  "edge_id",
  "selection_frequency",
  "selection_frequency_mc_se",
  "selection_frequency_lower_95",
  "selection_frequency_upper_95",
  "sign_stability",
  "sign_agreement_fraction",
  "cv_rsq",
  "active_edge"
), colnames(step1_sensitive$grn_result$stability_diagnostics)), drop = FALSE]
```

Targets with inadequate bootstrap completion or non-positive out-of-fold R-squared must not be interpreted as stable regulatory edges.

## 4. Structural-domain sensitivity: GREAT versus Signac

Changing the canonical GREAT domains to the narrower Signac-style domain rule
changes the biological candidate hypothesis and therefore the Pando structural
fingerprint:

```r
step1_signac <- rc_regcompass_step_grn(
  object = A,
  gem = gem,
  outdir = "RegCompass_restart/01_grn_signac",
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
  )
)
```

Both Pando methods are structural domain rules in this workflow; neither uses a
fitted target-expression correlation threshold to admit candidates. GREAT is
canonical because it retains distal regulatory hypotheses, while Signac remains
a narrower structural sensitivity analysis. Compare complete-GPR cores only
after acknowledging that differences can arise from the changed structural
candidate universe, not only from coefficient estimation.

## 5. Rebuild dependent Stage 3–6 objects

```r
step3_sensitive <- rc_regcompass_step_meta_modules(
  grn = step1_sensitive,
  metacells = step2,
  gem = gem,
  outdir = "RegCompass_restart/03_meta_modules_sensitive"
)

step4_sensitive <- rc_regcompass_step_layer1(
  metacells = step2,
  meta_modules = step3_sensitive,
  gem = gem,
  outdir = "RegCompass_restart/04_layer1_sensitive",
  regulatory_alpha = 1,
  gpr_and_method = "min",
  parallel = TRUE,
  BPPARAM = upstream_bp
)

step5_sensitive <- rc_regcompass_step_layer2(
  layer1 = step4_sensitive,
  meta_modules = step3_sensitive,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "RegCompass_restart/05_layer2_sensitive",
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

result_sensitive <- rc_regcompass_step_results(
  grn = step1_sensitive,
  metacells = step2,
  meta_modules = step3_sensitive,
  layer1 = step4_sensitive,
  layer2 = step5_sensitive,
  gem = gem,
  outdir = "RegCompass_restart/06_results_sensitive",
  species = "human"
)
```

Do not pair `step1_sensitive` with the old `step3`, because target genes and core reactions may have changed.

## 6. GPR aggregation sensitivity

```r
step4_mean <- rc_regcompass_step_layer1(
  metacells = step2,
  meta_modules = step3,
  gem = gem,
  outdir = "RegCompass_restart/04_layer1_mean",
  regulatory_alpha = 1,
  gpr_and_method = "mean",
  parallel = TRUE,
  BPPARAM = upstream_bp
)
```

The canonical default `min` represents limiting-subunit logic. `median` and `mean` are sensitivity analyses, not equivalent biological assumptions. After changing Stage 4, rebuild Stage 5 and Stage 6.

## 7. Medium sensitivity

```r
low_glucose_medium <- rc_make_medium_scenarios(
  gem,
  scenario = "low_glucose",
  species = "human"
)

step5_low_glucose <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = low_glucose_medium,
  outdir = "RegCompass_restart/05_layer2_low_glucose",
  model_mode = "meta_module_gem",
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    model_params = list(
      completion_time_limit = 900,
      fastcore_epsilon = 1e-4,
      max_support_reactions = 2500,
      strict = TRUE
    )
  ),
  parallel = TRUE,
  BPPARAM = layer2_bp
)
```

A new medium may require a different global FASTCORE support set. Structural comparability applies within a medium, not across different media.

```r
step5_low_glucose$model_cache_summary
```

`completion_time_limit` constrains union-GEM construction only. Directional scoring LPs do not accept a scoring time limit.

## 8. Handoff

Use the unchanged `step3`, `step4`, and `step5` for [Tutorial 4](tutorial-04-targeted-reaction-remapping.md).

Use either `result` or `result_sensitive` in [Tutorial 5](tutorial-05-condition-differential-analysis.md), but keep each result paired with its own stage checkpoints when tracing evidence.

See [Pando and multitask GRN parameter policy](grn-parameter-policy.md) for the canonical defaults and the boundary between structural and modelling sensitivity analyses.
