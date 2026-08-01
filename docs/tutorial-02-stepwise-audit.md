# Tutorial 2: stepwise workflow

Use the stepwise API to save, inspect and restart the canonical stages.
Equations are in [Tutorial 3](tutorial-03-mathematical-model.md). Public API:
[functions.md](functions.md).

## Detailed progress and audit log

All public stages accept `progress = TRUE`. Stage 1 reports input validation, design resolution, normalization, Pando runtime checks, target selection, ATAC filtering, candidate initialization, motif mapping, nested CV, contract extraction and artifact writing. The long-running nested-CV event includes the number of cell types, conditions, metabolic targets, outer/inner folds and lambda values. The RegCompass route skips exact motif-hit coordinates and retains only Pando's binary peak-by-motif incidence matrix, which is the representation used for candidate construction.

```r
options(RegCompassR.progress = TRUE)
step1 <- rc_regcompass_step_grn(..., progress = TRUE)
progress_log <- read.delim(
  "RegCompass_steps/01_grn/step_progress.tsv",
  check.names = FALSE
)
progress_log[, c("phase", "elapsed_hms", "detail", "context")]
```

`step_progress.tsv` is written even when `progress = FALSE`; that setting only suppresses console messages. Pando target-level messages are enabled automatically when Stage 1 progress is enabled. Errors terminate immediately and the last audit row is `stage_error`.

Both Stage 1 routes are instrumented. Condition-aware runs report the fused
nested-CV phases, while standard Pando runs print each cell-type job's candidate
initialization, binary motif mapping, GRN fit, edge extraction, and artifact
write. An error or user interrupt is printed immediately with its original
message and the active phase, even when routine progress output is disabled.

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

`BPPARAM = TRUE` is invalid. Do not set `parallel` inside
`pando_infer_args`.

## Stage 1: condition-aware or standard Pando

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
    min_cells = 100L,
    min_abs_estimate = 0,
    min_model_rsq = 0.1,
    pando_infer_args = list(
      candidate_screen = "motif_domain",
      condition_mix = 0.5,
      condition_weight = "equal",
      outer_nfolds = 5L,
      inner_nfolds = 5L,
      lambda_selection = "lambda.1se",
      scale = TRUE
    )
  ),
  parallel = TRUE,
  BPPARAM = upstream_bp
)
```

```r
step1$params$analysis_mode
step1$grn_result$condition_grn_fits
step1$grn_result$tf_peak_gene_condition
```

With fewer than two condition levels, Stage 1 uses `standard_pando` and No
condition coefficients are calculated.

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

Each cell type has one independent graph. All conditions are joint in that
graph, and final metacells are condition-pure. Set `overwrite = TRUE` after
changing cells, reductions, dimensions, gamma, seed or thresholds.

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

## Stage 4: condition-full reaction support

```r
step4 <- rc_regcompass_step_layer1(
  grn = step1,
  metacells = step2,
  meta_modules = step3,
  gem = gem,
  outdir = "RegCompass_steps/04_layer1",
  projection_component = "condition",
  comparison_support = "auto",
  regulatory_alpha = 1,
  gpr_and_method = "min",
  parallel = TRUE,
  BPPARAM = upstream_bp
)
```

```r
step4$gene_projection_condition_full_oof
step4$gene_projection_common_oof
step4$gene_projection_condition_unique_oof
step4$reaction_expression_condition_full_oof
step4$reaction_expression_common_oof
step4$projection_provenance
```

Condition-full OOF is primary. Common support is the jointly estimable
component. Each non-estimable edge side contributes zero.

## Build Stage 5 media

### Plasma scenarios

```r
human_medium <- rc_make_medium_scenarios(
  gem = human_gem,
  scenario = "normal_human_plasma",
  species = "human"
)

mouse_medium <- rc_make_medium_scenarios(
  gem = mouse_gem,
  scenario = "mouse_plasma",
  species = "mouse"
)
```

`normal_human_plasma` encodes HPLM composition from *Cell* 2017 and the updated
HPLM formulation from *Cell Metabolism* 2021. Plasmax from *Science Advances*
2019 is validation only and is not averaged with HPLM. `mouse_plasma` uses a
conservative metabolite set anchored to absolute plasma and interstitial-fluid
measurements in *Nature* 2026. Unsupported mouse components are omitted rather
than inferred from human HPLM.

### Culture challenge scenarios

```r
medium_scenarios <- rc_make_medium_scenarios(
  gem = human_gem,
  scenario = c(
    "high_glucose",
    "low_glucose",
    "high_lactate",
    "low_lactate",
    "low_glutamine"
  ),
  species = "human"
)
```

All five challenge scenarios use the identical HPLM 2017/2021 basal nutrient
composition. The challenge paper replaces only the named target concentration:

```text
high_glucose   glucose 25 mM    Han 2015
low_glucose    glucose 1 mM     Han 2015
high_lactate   lactate 20 mM    San-Millan 2020
low_lactate    lactate 0.5 mM   Cho 2025
low_glutamine  glutamine 0.5 mM Visagie 2015 Methods
```

This removes RPMI-versus-DMEM-versus-Plasmax composition as a confounder between
challenge scenarios. Inspect composition, validation and challenge provenance:

```r
unique(medium_scenarios[, intersect(c(
  "medium_scenario_id",
  "medium_background_id",
  "background_reference_label",
  "background_reference_doi",
  "background_validation_reference_label",
  "background_validation_reference_doi",
  "challenge_reference_label",
  "challenge_reference_doi",
  "scenario_construction"
), colnames(medium_scenarios))])
```

Target concentration-derived uptake caps remain sensitivity assumptions rather
than measured transporter rates.

### User-defined composition

```r
custom_medium <- data.frame(
  medium_scenario_id = "my_measured_medium",
  exchange_reaction_id = c("EX_glc_D_e", "EX_gln_L_e"),
  lb = c(-0.20, -0.10),
  ub = c(1, 1),
  available = TRUE,
  reference_label = "Optional experiment or publication label",
  reference_doi = "10.xxxx/optional.reference",
  stringsAsFactors = FALSE
)

custom_medium <- rc_make_medium_scenarios(
  gem = human_gem,
  scenario = "custom",
  species = "human",
  custom_medium = custom_medium
)
```

`scenario = NULL` is accepted for a custom-only run. `custom_metabolites` can be
used instead of reaction-level bounds. Built-in and custom scenarios may also be
returned together.

Changing the scenario list, custom composition, exchange limit, target
`uptake_scale`, or any resulting bound invalidates Stage 5 and downstream
results but does not require rerunning Stages 1-4. See
[Medium scenarios and evidence](medium-presets.md).

## Stage 5: shared model and LP scoring

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

```r
step5$penalty_condition_full_oof
step5$penalty_common_oof
step5$penalty_condition_unique_increment
step5$penalty_rna_only
step5$model_cache_summary
```

All routes use the exact same medium-specific model, bounds, target directions
and `vmax`. The five retired guardrails are absent from the result schema.

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
result$common_support_component_summary
result$condition_unique_penalty_increment_summary
```

## Optional targeted reaction remapping

```r
targeted <- rc_regcompass_step_target_union(
  layer1 = step4,
  meta_modules = step3,
  layer2 = step5,
  gem = gem,
  outdir = "RegCompass_steps/07_targeted_remapping",
  core_reaction_ids = c("MAR04381", "MAR04379"),
  layer2_args = list(target_direction = "both", solver = "highs"),
  parallel = TRUE,
  BPPARAM = layer2_bp
)
```

This pass uses `step4$reaction_expression`, which is the condition-full primary
reaction-expression route in condition mode, and reuses the exact cached Stage 5
union GEM without rebuilding it or rerunning FASTCORE.

Saved stages:

```text
01_grn/step_grn.rds
02_metacells/step_metacells.rds
03_meta_modules/step_meta_modules.rds
04_layer1/step_layer1.rds
05_layer2/step_layer2.rds
06_results/regcompass_result.rds
07_targeted_remapping/step_target_union.rds
```

See [Tutorial 4](tutorial-04-targeted-reaction-remapping.md) for targeted
reaction extension and [Tutorial 5](tutorial-05-condition-differential-analysis.md)
for condition statistics.
