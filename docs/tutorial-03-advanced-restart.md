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
- multitask penalties, CV rule, bootstrap number, or active-edge thresholds;
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

## 3. GRN bootstrap sensitivity

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
      peak_to_gene_method = "Signac",
      min_tf_detection = 0.02,
      min_peak_detection = 0.02,
      min_target_detection = 0.02
    )
  ),
  multitask_args = list(
    alpha = 0.25,
    global_penalty_factor = 1,
    deviation_penalty_factor = 3,
    lambda_rule = "lambda.1se",
    nfolds = 5,
    n_bootstrap = 200,
    min_selection_frequency = 0.8,
    min_sign_stability = 0.9,
    min_bootstrap_success_fraction = 0.8,
    candidate_screen_threshold = 0,
    max_edges_per_target = Inf,
    seed = 12345L
  )
)
```

`alpha` must remain strictly between zero and one. A lasso component is required for sparse bootstrap selection, and a ridge component is required for a unique symmetric condition-deviation solution.

Increasing `n_bootstrap` reduces Monte Carlo error; it does not add independent biological replicates.

Compare the same structural universe before comparing active-edge counts:

```r
old_universe <- unique(
  step1$grn_result$tf_peak_gene_candidates$edge_universe_id
)
new_universe <- unique(
  step1_sensitive$grn_result$tf_peak_gene_candidates$edge_universe_id
)

old_universe
new_universe
```

If structural thresholds changed, different universes are expected. Then changes in active-edge counts reflect both candidate-space and coefficient-selection changes.

Check bootstrap completion:

```r
with(
  step1_sensitive$grn_result$target_model_diagnostics,
  summary(n_bootstrap_success / n_bootstrap_requested)
)
```

Targets with inadequate bootstrap completion must not be interpreted as biologically unstable edges.

## 4. Rebuild dependent Stage 3–6 objects

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

## 5. GPR aggregation sensitivity

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

## 6. Medium sensitivity

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

## 7. Handoff

Use the unchanged `step3`, `step4`, and `step5` for [Tutorial 4](tutorial-04-targeted-reaction-remapping.md).

Use either `result` or `result_sensitive` in [Tutorial 5](tutorial-05-condition-differential-analysis.md), but keep each result paired with its own stage checkpoints when tracing evidence.
