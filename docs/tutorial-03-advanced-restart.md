# Tutorial Level 3: restart, sensitivity, and diagnostics

RegCompass saves classed stage objects so downstream stages can be rerun without repeating unchanged work. Incompatible workflow settings, GEM fingerprints, and ordered-unit contracts are rejected.

## Load a completed stepwise run

```r
step1 <- readRDS("RegCompass_steps/01_grn/step_grn.rds")
step2 <- readRDS("RegCompass_steps/02_metacells/step_metacells.rds")
step3 <- readRDS("RegCompass_steps/03_meta_modules/step_meta_modules.rds")
step4 <- readRDS("RegCompass_steps/04_layer1/step_layer1.rds")
step5 <- readRDS("RegCompass_steps/05_layer2/step_layer2.rds")
```

## Restart boundaries

### Rerun Stage 1 onward

Rerun Stage 1 and every downstream stage after changing:

- motif matrices, genome, or Pando regulatory regions;
- `pando_design_args`, including peak-to-gene domains or structural detection filters;
- `multitask_args`, including elastic-net penalties, lambda rule, condition-stratified folds, `n_bootstrap`, or active-edge thresholds;
- `grn_mode`;
- condition or cell-type metadata;
- single-cell RNA or ATAC matrices;
- GEM GPR target genes.

The canonical workflow does not accept a biological-sample column. Stage 1 centring, CV, and bootstrap are condition based.

Legacy `padj_threshold`, `min_abs_estimate`, `min_model_rsq`, and `pando_infer_args` apply only under:

```r
grn_mode = "legacy_condition_pando"
```

### Rerun Stage 2 onward

Rerun metacells and downstream stages after changing:

- cell membership;
- RNA/ATAC matrices used for aggregation;
- condition or cell-type metadata;
- RNA/ATAC reductions or dimensions;
- `gamma`, seed, or metacell thresholds.

Condition remains the only hard stratum. Cell type is the SuperCell2 label. No sample balancing or sample composition is calculated.

Changing Stage 1 without changing cells does not require Stage 2 to be rebuilt, but Stage 3 onward must be rerun because condition target genes may change.

### Rerun Stage 3 onward

Rerun condition meta-modules and downstream stages after changing:

- active bootstrap-stable condition edges or target genes;
- a custom `meta_module_args$subsystem_table`;
- subsystem, KEGG, Reactome, or master-Rhea annotations;
- GEM GPR rules.

Stage 3 always performs one ordered pass:

```text
complete-GPR cores
→ core subsystems
→ direct KEGG/Reactome equivalence
→ direct master-Rhea equivalence
→ stop
```

### Rerun Stage 4 onward

Rerun Layer 1 and downstream stages after changing:

- `regulatory_alpha`;
- `gpr_and_method` among `min`, `median`, and `mean`;
- gene or ATAC half-saturation options;
- metacell RNA or ATAC evidence;
- Stage 1 stable coefficients or ATAC projection weights.

The default GPR-AND method is `min`.

### Rerun Stage 5 onward

Rerun Layer 2 and results after changing:

- medium composition or exchange bounds;
- `target_direction`, `omega`, solver, or `flux_threshold`;
- `completion_time_limit`, `fastcore_epsilon`, `max_support_reactions`, or `strict`;
- Stage 3 merged reaction membership;
- Stage 4 reaction penalties.

`completion_time_limit` applies only to construction of the medium-specific union GEM. Scoring LPs do not accept a time-limit parameter.

## Stage 1 bootstrap sensitivity example

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
    candidate_screen_threshold = 0,
    max_edges_per_target = Inf,
    seed = 12345L
  )
)
```

Increasing `n_bootstrap` reduces Monte Carlo error in `selection_frequency` and `sign_stability`; it does not add independent biological replicates. Compare candidate and active-edge counts only after confirming shared `edge_universe_id` values:

```r
step1$grn_result$celltype_fit_status
step1_sensitive$grn_result$celltype_fit_status

unique(step1$grn_result$tf_peak_gene_candidates$edge_universe_id)
unique(step1_sensitive$grn_result$tf_peak_gene_candidates$edge_universe_id)

table(step1$grn_result$tf_peak_gene_significant$Group)
table(step1_sensitive$grn_result$tf_peak_gene_significant$Group)
```

Inspect bootstrap completion before interpreting stability:

```r
with(
  step1_sensitive$grn_result$target_model_diagnostics,
  summary(n_bootstrap_success / n_bootstrap_requested)
)
```

Targets with failed bootstrap fits should be diagnosed rather than treated as low biological stability.

## Rerun Stage 4 with another GPR-AND rule

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

Use `median` or `mean` for sensitivity analysis. The canonical default remains `min`.

## Rebuild Stage 5 under another medium

```r
new_medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "low_glucose",
  species = "human"
)

step5_new <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = new_medium_scenarios,
  outdir = "RegCompass_restart/05_layer2",
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

The new medium may require a different global FASTCORE support set. Within each medium, all conditions still share one final model.

## Diagnose model completion

```r
summary <- step5_new$model_cache_summary
summary[, c(
  "medium_scenario",
  "n_merged_biological_reactions",
  "n_global_fastcore_support_reactions",
  "target_status",
  "file",
  "file_checksum"
)]

model <- readRDS(summary$file[[1]])
table(model$closure_diagnostics$completion_status)
```

Common statuses include `already_feasible`, `global_fastcore_completed`, `parent_blocked`, `unresolved`, and `no_allowed_direction`.
